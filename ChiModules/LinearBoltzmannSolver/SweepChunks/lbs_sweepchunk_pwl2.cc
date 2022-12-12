#include "lbs_sweepchunk_pwl.h"

#include "chi_runtime.h"
#include "LinearBoltzmannSolver/Groupset/lbs_groupset.h"

#ifdef WITH_READABLE_CHUNK

double* lbs::SweepChunkPWL::Upwinder::
  GetUpwindPsi(int fj, bool local, bool boundary) const
{
  double* psi;
  if (local)             psi = fluds->UpwindPsi(spls_index,
                                                in_face_counter,
                                                fj,0,angle_set_index);
  else if (not boundary) psi = fluds->NLUpwindPsi(preloc_face_counter,
                                                  fj,0,angle_set_index);
  else                   psi = angle_set->PsiBndry(bndry_id,
                                                   angle_num,
                                                   cell_local_id,
                                                   f, fj, gs_gi, gs_ss_begin,
                                                   surface_source_active);
  return psi;
}

double* lbs::SweepChunkPWL::Upwinder::
  GetDownwindPsi(int fi, bool local, bool boundary, bool reflecting_bndry) const
{
  double* psi;
  if (local)                 psi = fluds->
    OutgoingPsi(spls_index,
                out_face_counter,
                fi, angle_set_index);
  else if (not boundary)     psi = fluds->
    NLOutgoingPsi(deploc_face_counter,
                  fi, angle_set_index);
  else if (reflecting_bndry) psi = angle_set->
    ReflectingPsiOutBoundBndry(bndry_id, angle_num,
                               cell_local_id, f,
                               fi, gs_ss_begin);
  else
    psi = nullptr;

  return psi;
}

//###################################################################
/**Actual sweep function*/
void lbs::SweepChunkPWL::
  Sweep(chi_mesh::sweep_management::AngleSet *angle_set)
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
  std::vector<double>& output_phi = GetDestinationPhi();
  std::vector<double>& output_psi = GetDestinationPsi();

  const SubSetInfo& grp_ss_info = groupset.grp_subset_infos[angle_set->ref_subset];
  const size_t gs_ss_size  = grp_ss_info.ss_size;
  const size_t gs_ss_begin = grp_ss_info.ss_begin;
  const int    gs_gi = groupset.groups[gs_ss_begin].id; // Groupset subset first group number

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto const& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();
  auto const& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();

  // ========================================================== Loop over each cell
  size_t num_loc_cells = spds->spls.item_id.size();
  for (size_t spls_index = 0; spls_index < num_loc_cells; ++spls_index)
  {
    const int cell_local_id = spds->spls.item_id[spls_index];
    const auto& cell = grid_view->local_cells[cell_local_id];
    const auto& fe_intgrl_values = grid_fe_view.GetUnitIntegrals(cell);
    const auto num_faces = cell.faces.size();
    const int num_nodes = static_cast<int>(fe_intgrl_values.NumNodes());
    auto& transport_view = grid_transport_view[cell.local_id];
    const auto& sigma_tg = transport_view.XS().sigma_t;
    std::vector<bool> face_incident_flags(num_faces, false);
    std::vector<double> face_mu_values(num_faces, 0.0);

    // =================================================== Get Cell matrices
    const auto& G           = fe_intgrl_values.GetIntV_shapeI_gradshapeJ();
    const auto& M           = fe_intgrl_values.GetIntV_shapeI_shapeJ();
    const auto& M_surf      = fe_intgrl_values.GetIntS_shapeI_shapeJ();
    const auto& IntS_shapeI = fe_intgrl_values.GetIntS_shapeI();

    // =================================================== Loop over angles in set
    const int ni_deploc_face_counter = deploc_face_counter;
    const int ni_preloc_face_counter = preloc_face_counter;
    const size_t as_num_angles = angle_set->angles.size();
    for (size_t angle_set_index = 0; angle_set_index<as_num_angles; ++angle_set_index)
    {
      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;
      const int angle_num = angle_set->angles[angle_set_index];
      const chi_mesh::Vector3& omega = groupset.quadrature->omegas[angle_num];
      const double wt = groupset.quadrature->weights[angle_num];

      // ============================================ Gradient matrix
      for (int i = 0; i < num_nodes; ++i)
        for (int j = 0; j < num_nodes; ++j)
          Amat[i][j] = omega.Dot(G[i][j]);

      for (int gsg = 0; gsg < gs_ss_size; ++gsg)
        b[gsg].assign(num_nodes, 0.0);

      // ============================================ Upwinding structure
      Upwinder upwind{fluds, angle_set, spls_index, angle_set_index,
                      /*in_face_counter*/0,
                      /*preloc_face_counter*/0,
                      /*out_face_counter*/0,
                      /*deploc_face_counter*/0,
                      /*bndry_id*/0, angle_num, cell.local_id,
                      /*f*/0,gs_gi, gs_ss_begin,
                      surface_source_active};

      // ============================================ Surface integrals
      int in_face_counter = -1;
      for (int f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces[f];
        const double mu = omega.Dot(face.normal);
        face_mu_values[f] = mu;

        if (mu >= 0.0) continue;

        face_incident_flags[f] = true;
        const bool local = transport_view.IsFaceLocal(f);
        const bool boundary = not face.has_neighbor;
        const uint64_t bndry_id = face.neighbor_id;

        if (local)             ++in_face_counter;
        else if (not boundary) ++preloc_face_counter;

        upwind.in_face_counter = in_face_counter;
        upwind.preloc_face_counter = preloc_face_counter;
        upwind.bndry_id = bndry_id;
        upwind.f = f;

        const size_t num_face_indices = face.vertex_ids.size();
        for (int fi = 0; fi < num_face_indices; ++fi)
        {
          const int i = fe_intgrl_values.FaceDofMapping(f,fi);
          for (int fj = 0; fj < num_face_indices; ++fj)
          {
            const int j = fe_intgrl_values.FaceDofMapping(f,fj);

            const double* psi = upwind.GetUpwindPsi(fj, local, boundary);

            const double mu_Nij = -mu * M_surf[f][i][j];
            Amat[i][j] += mu_Nij;
            for (int gsg = 0; gsg < gs_ss_size; ++gsg)
              b[gsg][i] += psi[gsg]*mu_Nij;

          }//for face j
        }//for face i
      } // for f

      // ========================================== Looping over groups
      for (int gsg = 0; gsg < gs_ss_size; ++gsg)
      {
        const int g = gs_gi+gsg;

        // ============================= Contribute source moments
        // q = M_n^T * q_moms
        for (int i = 0; i < num_nodes; ++i)
        {
          double temp_src = 0.0;
          for (int m = 0; m < num_moms; ++m)
          {
            const size_t ir = transport_view.MapDOF(i, m, g);
            temp_src += m2d_op[m][angle_num]*q_moments[ir];
          }//for m
          source[i] = temp_src;
        }//for i

        // ============================= Mass Matrix and Source
        // Atemp  = Amat + sigma_tgr * M
        // b     += M * q
        const double sigma_tgr = sigma_tg[g];
        for (int i = 0; i < num_nodes; ++i)
        {
          double temp = 0.0;
          for (int j = 0; j < num_nodes; ++j)
          {
            const double Mij = M[i][j];
            Atemp[i][j] = Amat[i][j] + Mij*sigma_tgr;
            temp += Mij*source[j];
          }//for j
          b[gsg][i] += temp;
        }//for i

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
            output_phi[ir + gsg] += wn_d2m * b[gsg][i];
        }
      }

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
            output_psi[ir + gsg] = b[gsg][i];
        }//for i
      }//if save psi

      int out_face_counter = -1;
      for (int f = 0; f < num_faces; ++f)
      {
        if (face_incident_flags[f]) continue;
        double mu = face_mu_values[f];

        // ============================= Set flags and counters
        out_face_counter++;
        const auto& face = cell.faces[f];
        const bool local = transport_view.IsFaceLocal(f);
        const bool boundary = not face.has_neighbor;
        const size_t num_face_indices = face.vertex_ids.size();
        const std::vector<double>& IntF_shapeI = IntS_shapeI[f];
        const uint64_t bndry_id = face.neighbor_id;

        bool reflecting_bndry = false;
        if (boundary)
          if (angle_set->ref_boundaries[bndry_id]->IsReflecting())
            reflecting_bndry = true;

        if (not boundary and not local) ++deploc_face_counter;

        upwind.out_face_counter = out_face_counter;
        upwind.deploc_face_counter = deploc_face_counter;
        upwind.bndry_id = bndry_id;
        upwind.f = f;

        for (int fi = 0; fi < num_face_indices; ++fi)
        {
          const int i = fe_intgrl_values.FaceDofMapping(f,fi);

          double* psi = upwind.GetDownwindPsi(fi, local, boundary, reflecting_bndry);

          if (not boundary or reflecting_bndry)
            for (int gsg = 0; gsg < gs_ss_size; ++gsg)
              psi[gsg] = b[gsg][i];
          if (boundary and not reflecting_bndry)
            for (int gsg = 0; gsg < gs_ss_size; ++gsg)
              transport_view.AddOutflow(gs_gi + gsg,
                                        wt*mu*b[gsg][i]*IntF_shapeI[i]);
        }//for fi
      }//for face
    } // for n
  } // for cell
}//Sweep

#endif