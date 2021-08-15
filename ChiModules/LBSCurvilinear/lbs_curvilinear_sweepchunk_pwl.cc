#include "LBSCurvilinear/lbs_curvilinear_sweepchunk_pwl.h"

#include "ChiMath/Quadratures/curvilinear_angular_quadrature.h"


LBSCurvilinear::SweepChunkPWL::
  SweepChunkPWL(std::shared_ptr<chi_mesh::MeshContinuum> grid_ptr,
                SpatialDiscretization_PWLD& discretization_primary,
                SpatialDiscretization_PWLD& discretization_secondary,
                std::vector<LinearBoltzmann::CellLBSView>& cell_transport_views,
                std::vector<double>& destination_phi,
                std::vector<double>& destination_psi,
                const std::vector<double>& source_moments,
                LBSGroupset& in_groupset,
                const TCrossSections& in_xsections,
                const int in_num_moms,
                const int in_max_num_cell_dofs)
  : LinearBoltzmann::SweepChunkPWL(
                     std::move(grid_ptr),
                     discretization_primary,
                     cell_transport_views,
                     destination_phi,
                     destination_psi,
                     source_moments,
                     in_groupset,
                     in_xsections,
                     in_num_moms,
                     in_max_num_cell_dofs)
  , grid_fe_view_secondary(discretization_secondary)
  , unknown_manager()
  , psi_start()
  , psi_sweep()
  , normal_vector_boundary()
{
  const auto curvilinear_product_quadrature =
    std::dynamic_pointer_cast<chi_math::CurvilinearAngularQuadrature>(groupset.quadrature);

  if (!curvilinear_product_quadrature)
    throw std::invalid_argument("LBSCurvilinear::SweepChunkPWL::SweepChunkPWL : "
                                "invalid angular quadrature");

  //  configure unknown manager for quantities that depend on polar level
  for (const auto& dir_set : curvilinear_product_quadrature->GetDirectionMap())
    unknown_manager.AddUnknown(chi_math::UnknownType::VECTOR_N,
                               groupset.groups.size());

  //  allocate storage for starting direction and for sweeping dependency
  const unsigned int n_dof =
    discretization_primary.GetNumLocalDOFs(unknown_manager);
  psi_start.resize(n_dof);
  psi_sweep.resize(n_dof);

  //  initialise mappings from direction linear index
  for (const auto& dir_set : curvilinear_product_quadrature->GetDirectionMap())
    for (const auto& dir_idx : dir_set.second)
    {
      map_polar_level.emplace(dir_idx, dir_set.first);

      const auto dir_start = (dir_idx == dir_set.second.front());
      const auto dir_final = (dir_idx == dir_set.second.back());
      map_start_final.emplace(dir_idx, std::make_pair(dir_start, dir_final));
    }

  //  set normal vector for symmetric boundary condition
  const int d =
    (grid_view->local_cells[0].Type() == chi_mesh::CellType::SLAB) ? 2 : 0;
  normal_vector_boundary = chi_mesh::Vector3(0.0, 0.0, 0.0);
  normal_vector_boundary(d) = 1;
}


void
LBSCurvilinear::SweepChunkPWL::Sweep(chi_mesh::sweep_management::AngleSet* angle_set)
{
  if (!a_and_b_initialized)
  {
    Amat.resize(max_num_cell_dofs, std::vector<double>(max_num_cell_dofs));
    Atemp.resize(max_num_cell_dofs, std::vector<double>(max_num_cell_dofs));
    b.resize(num_grps, std::vector<double>(max_num_cell_dofs, 0.0));
    source.resize(max_num_cell_dofs, 0.0);
    a_and_b_initialized = true;
  }

  const auto spds = angle_set->GetSPDS();
  const auto fluds = angle_set->fluds;
  const bool surface_source_active = IsSurfaceSourceActive();
  std::vector<double>& output_vector = GetDestinationPhi();

  const GsSubSet& subset = groupset.grp_subsets[angle_set->ref_subset];
  const int gs_ss_size  = groupset.grp_subset_sizes[angle_set->ref_subset];
  const int gs_ss_begin = subset.first;
  const int gs_gi = groupset.groups[gs_ss_begin].id; // Groupset subset first group number

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  const auto& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();
  const auto& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();

  const auto curvilinear_product_quadrature =
    std::dynamic_pointer_cast<chi_math::CurvilinearAngularQuadrature>(groupset.quadrature);

  //========================================================== Loop over each cell
  size_t num_loc_cells = spds->spls.item_id.size();
  for (size_t spls_index = 0; spls_index < num_loc_cells; ++spls_index)
  {
    const int cell_local_id = spds->spls.item_id[spls_index];
    const auto& cell = grid_view->local_cells[cell_local_id];
    const auto& fe_intgrl_values = grid_fe_view.GetUnitIntegrals(cell);
    const auto& fe_intgrl_values_secondary = grid_fe_view_secondary.GetUnitIntegrals(cell);
    const auto num_faces = cell.faces.size();
    const int num_nodes = static_cast<int>(fe_intgrl_values.NumNodes());
    auto& transport_view = grid_transport_view[cell.local_id];
    const int xs_mapping = transport_view.XSMapping();
    const auto& sigma_tg = xsections[xs_mapping]->sigma_t;
    std::vector<bool> face_incident_flags(num_faces, false);
    std::vector<double> face_mu_values(num_faces, 0.0);

    //=================================================== Get Cell matrices
    const auto& G      = fe_intgrl_values.GetIntV_shapeI_gradshapeJ();
    const auto& M      = fe_intgrl_values.GetIntV_shapeI_shapeJ();
    const auto& M_surf = fe_intgrl_values.GetIntS_shapeI_shapeJ();

    const auto& Maux   = fe_intgrl_values_secondary.GetIntV_shapeI_shapeJ();


    //=================================================== Loop over angles in set
    const int ni_deploc_face_counter = deploc_face_counter;
    const int ni_preloc_face_counter = preloc_face_counter;
    const size_t as_num_angles = angle_set->angles.size();
    for (size_t angle_set_index = 0; angle_set_index < as_num_angles; ++angle_set_index)
    {
      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;
      const auto& angle_num = angle_set->angles[angle_set_index];
      const auto& omega = groupset.quadrature->omegas[angle_num];

      const auto polar_level = map_polar_level[angle_num];
      const auto start_direction = map_start_final[angle_num].first;
    //const auto final_direction = map_start_final[angle_num].second;

      const auto& fac_diamond_difference =
        curvilinear_product_quadrature->GetDiamondDifferenceFactor()[angle_num];
      const auto& fac_streaming_operator =
        curvilinear_product_quadrature->GetStreamingOperatorFactor()[angle_num];


      // ============================================ Gradient matrix
      for (size_t i = 0; i < num_nodes; ++i)
        for (size_t j = 0; j < num_nodes; ++j)
          Amat[i][j] = omega.Dot(G[i][j]) + fac_streaming_operator * Maux[i][j];


      // ============================================ Source initialization
      for (int gsg = 0; gsg < gs_ss_size; ++gsg)
        b[gsg].assign(num_nodes, 0);

      for (size_t i = 0; i < num_nodes; ++i)
        for (size_t j = 0; j < num_nodes; ++j)
        {
          const auto jr = grid_fe_view.MapDOFLocal(cell, j, unknown_manager, polar_level, gs_gi);
          for (int gsg = 0; gsg < gs_ss_size; ++gsg)
            b[gsg][i] += fac_streaming_operator * Maux[i][j] * psi_sweep[jr+gsg];
        }


      // ============================================ Surface integrals
      int in_face_counter = -1;
      for (int f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces[f];
        const double mu = omega.Dot(face.normal);

        if (mu < 0.0) // Upwind
        {
          face_incident_flags[f] = true;
          const bool local = transport_view.IsFaceLocal(f);
          const bool boundary = not face.has_neighbor;
          const size_t num_face_indices = face.vertex_ids.size();
          if (local)
          {
            in_face_counter++;
            for (int fi = 0; fi < num_face_indices; ++fi)
            {
              const int i = fe_intgrl_values.FaceDofMapping(f,fi);
              for (int fj = 0; fj < num_face_indices; ++fj)
              {
                const int j = fe_intgrl_values.FaceDofMapping(f,fj);
                const double *psi = fluds->UpwindPsi(spls_index,in_face_counter,fj,0,angle_set_index);
                const double mu_Nij = -mu * M_surf[f][i][j];
                Amat[i][j] += mu_Nij;
                for (int gsg = 0; gsg < gs_ss_size; ++gsg)
                  b[gsg][i] += psi[gsg]*mu_Nij;
              }
            }
          }
          else if (not boundary)
          {
            preloc_face_counter++;
            for (int fi = 0; fi < num_face_indices; ++fi)
            {
              const int i = fe_intgrl_values.FaceDofMapping(f,fi);
              for (int fj = 0; fj < num_face_indices; ++fj)
              {
                const int j = fe_intgrl_values.FaceDofMapping(f,fj);
                const double *psi = fluds->NLUpwindPsi(preloc_face_counter,fj,0,angle_set_index);
                const double mu_Nij = -mu * M_surf[f][i][j];
                Amat[i][j] += mu_Nij;
                for (int gsg = 0; gsg < gs_ss_size; ++gsg)
                  b[gsg][i] += psi[gsg]*mu_Nij;
              }
            }
          }
          else
          {
            //------------------------------------------------------------------
            //  determine whether incoming direction is incident on the point
            //  of symmetry or on the axis of symmetry
            //  N.B.: a face is considered to be on the point/axis of symmetry
            //  if all are true:
            //    1. the face normal is antiparallel to $\vec{e}_{d}$
            //    2. all vertices of the face exhibit $v_{d} = 0$
            //  with $d = 2$ for 1D geometries and $d = 0$ for 2D geometries.
            //  Thanks to the verifications performed during initialisation,
            //  at this point it is necessary to confirm only the orientation.
            const bool incident_on_symmetric_boundary =
              (face.normal.Dot(normal_vector_boundary) < -0.999999);
            if (!incident_on_symmetric_boundary)
            {
              const uint64_t bndry_index = face.neighbor_id;
              for (int fi = 0; fi < num_face_indices; ++fi)
              {
                const int i = fe_intgrl_values.FaceDofMapping(f,fi);
                for (int fj = 0; fj < num_face_indices; ++fj)
                {
                  const int j = fe_intgrl_values.FaceDofMapping(f,fj);
                  const double *psi = angle_set->PsiBndry(bndry_index,
                                                          angle_num,
                                                          cell.local_id,
                                                          f, fj, gs_gi, gs_ss_begin,
                                                          surface_source_active);
                  const double mu_Nij = -mu * M_surf[f][i][j];
                  Amat[i][j] += mu_Nij;
                  for (int gsg = 0; gsg < gs_ss_size; ++gsg)
                    b[gsg][i] += psi[gsg]*mu_Nij;
                }
              }
            }
            else
            {
              for (int fi = 0; fi < num_face_indices; ++fi)
              {
                const int i = fe_intgrl_values.FaceDofMapping(f,fi);
                for (int fj = 0; fj < num_face_indices; ++fj)
                {
                  const int j = fe_intgrl_values.FaceDofMapping(f,fj);
                  const auto jr = grid_fe_view.MapDOFLocal(cell, j, unknown_manager, polar_level, gs_gi);
                  const double* psi = &psi_start[jr];
                  const double mu_Nij = -mu * M_surf[f][i][j];
                  Amat[i][j] += mu_Nij;
                  for (int gsg = 0; gsg < gs_ss_size; ++gsg)
                    b[gsg][i] += psi[gsg]*mu_Nij;
                }
              }
            }
            //------------------------------------------------------------------
          }
        } // if upwind
      } // for f



      // ========================================== Looping over groups
      for (int gsg = 0; gsg < gs_ss_size; ++gsg)
      {
        const auto g = gs_gi+gsg;

        // ============================= Contribute source moments
        for (int i = 0; i < num_nodes; ++i)
        {
          source[i] = 0;
          for (int m = 0; m < num_moms; ++m)
          {
            const size_t ir = transport_view.MapDOF(i, m, g);
            source[i] += m2d_op[m][angle_num] * q_moments[ir];
          }
        }

        // ============================= Mass Matrix and Source
        const auto& sigma_tgr = sigma_tg[g];
        for (size_t i = 0; i < num_nodes; ++i)
          for (size_t j = 0; j < num_nodes; ++j)
          {
            Atemp[i][j] = Amat[i][j] + M[i][j] * sigma_tgr;
            b[gsg][i] += M[i][j] * source[j];
          }

        // ============================= Solve system
        chi_math::GaussElimination(Atemp, b[gsg], num_nodes);
      }


      // ============================= Accumulate flux
      for (int m = 0; m < num_moms; ++m)
      {
        const double wn_d2m = d2m_op[m][angle_num];
        for (int i = 0; i < num_nodes; ++i)
        {
          const size_t ir = transport_view.MapDOF(i, m, gs_gi);
          for (int gsg = 0; gsg < gs_ss_size; ++gsg)
            output_vector[ir + gsg] += wn_d2m*b[gsg][i];
        }
      }


      // ============================= Moment callbacks
      for (auto& callback : moment_callbacks)
        callback(this, angle_set);


      // ============================= Save angular fluxes if needed
      if (save_angular_flux)
      {
        const auto& psi_uk_man = groupset.psi_uk_man;
        for (int i = 0; i < num_nodes; ++i)
        {
          int64_t ir = grid_fe_view.MapDOFLocal(cell,i,psi_uk_man,angle_num,0);
          for (int gsg = 0; gsg < gs_ss_size; ++gsg)
            psi_new_local[ir + gsg] = b[gsg][i];
        }
      }//if save psi


      //============================================= Outgoing fluxes
      int out_face_counter = -1;
      for (int f = 0; f < num_faces; ++f)
      {
        if (face_incident_flags[f]) continue;

        // ============================= Set flags and counters
        out_face_counter++;
        const auto& face = cell.faces[f];
        const bool local = transport_view.IsFaceLocal(f);
        const bool boundary = not face.has_neighbor;
        const size_t num_face_indices = face.vertex_ids.size();

        if (local)
        {
          for (int fi = 0; fi < num_face_indices; ++fi)
          {
            const int i = fe_intgrl_values.FaceDofMapping(f,fi);
            double *psi = fluds->OutgoingPsi(spls_index, out_face_counter, fi, angle_set_index);
            for (int gsg = 0; gsg < gs_ss_size; ++gsg)
              psi[gsg] = b[gsg][i];
          }
        }
        else if (not boundary)
        {
          deploc_face_counter++;
          for (int fi = 0; fi < num_face_indices; ++fi)
          {
            const int i = fe_intgrl_values.FaceDofMapping(f,fi);
            double *psi = fluds->NLOutgoingPsi(deploc_face_counter, fi, angle_set_index);
            for (int gsg = 0; gsg < gs_ss_size; ++gsg)
              psi[gsg] = b[gsg][i];
          }
        }
        else // Store outgoing reflecting Psi
        {
          const uint64_t bndry_index = face.neighbor_id;
          if (angle_set->ref_boundaries[bndry_index]->IsReflecting())
          {
            for (int fi = 0; fi < num_face_indices; ++fi)
            {
              const int i = fe_intgrl_values.FaceDofMapping(f,fi);
              double *psi = angle_set->ReflectingPsiOutBoundBndry(bndry_index, angle_num,
                                                                  cell.local_id, f,
                                                                  fi, gs_ss_begin);
              for (int gsg = 0; gsg < gs_ss_size; ++gsg)
                psi[gsg] = b[gsg][i];
            }
          }
        }
      }


      //  store starting direction angular intensity for each polar level
      if (start_direction)
        for (size_t i = 0; i < num_nodes; ++i)
        {
          const auto ir = grid_fe_view.MapDOFLocal(cell, i, unknown_manager, polar_level, gs_gi);
          for (int gsg = 0; gsg < gs_ss_size; ++gsg)
            psi_start[ir+gsg] = b[gsg][i];
        }

      //  update sweeping dependency angular intensity for each polar level
      //  (incoming for next interval)
      const auto f0 = 1/fac_diamond_difference;
      const auto f1 = f0 - 1;
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const auto ir = grid_fe_view.MapDOFLocal(cell, i, unknown_manager, polar_level, gs_gi);
        for (int gsg = 0; gsg < gs_ss_size; ++gsg)
          psi_sweep[ir+gsg] = f0 * b[gsg][i] - f1 * psi_sweep[ir+gsg];
      }


    }//for n

  }// for cell

}
