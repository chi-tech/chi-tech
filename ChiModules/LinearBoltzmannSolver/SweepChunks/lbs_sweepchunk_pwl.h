#ifndef LBS_SWEEPCHUNK_PWL_H
#define LBS_SWEEPCHUNK_PWL_H

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_polyhedron.h"

#include "ChiMesh/Cell/cell.h"
#include "ChiPhysics/chi_physics.h"

#include "ChiMath/chi_math.h"
#include "ChiModules/LinearBoltzmannSolver/GroupSet/lbs_groupset.h"
#include "ChiModules/LinearBoltzmannSolver/lbs_linear_boltzmann_solver.h"
#include "ChiMath/Quadratures/product_quadrature.h"

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiMesh/SweepUtilities/AngleAggregation/angleaggregation.h"

#include "ChiTimer/chi_timer.h"

#include "chi_mpi.h"
#include "chi_log.h"

extern ChiMath& chi_math_handler;
extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

typedef std::vector<std::shared_ptr<chi_physics::TransportCrossSections>> TCrossSections;

// ###################################################################
// Sweep chunk to compute the fixed source.
class LBSSweepChunkPWL : public chi_mesh::sweep_management::SweepChunk
{
protected:
  const std::shared_ptr<chi_mesh::MeshContinuum> grid_view;
  SpatialDiscretization_PWLD& grid_fe_view;
  const std::vector<LinearBoltzmann::CellLBSView>& grid_transport_view;
  const std::vector<double>* q_moments;
        LBSGroupset& groupset;
  const TCrossSections& xsections;
  const int num_moms;
  const int G;
  const int max_num_cell_dofs;

  //Runtime params
  bool a_and_b_initialized;
  std::vector<std::vector<double>> Amat;
  std::vector<std::vector<double>> Atemp;
  std::vector<double> source;

public:
  std::vector<std::vector<double>> b;
  int spls_index, angle_set_index;

  // ################################################## Constructor
  LBSSweepChunkPWL(std::shared_ptr<chi_mesh::MeshContinuum> grid_ptr,
                   SpatialDiscretization_PWLD& discretization,
                   const std::vector<LinearBoltzmann::CellLBSView>& cell_transport_views,
                   std::vector<double>* destination_phi,
                   const std::vector<double>* source_moments,
                         LBSGroupset& in_groupset,
                   const TCrossSections& in_xsections,
                   const int in_num_moms,
                   const int in_max_num_cell_dofs)
    : SweepChunk(destination_phi, false),
      grid_view(std::move(grid_ptr)),
      grid_fe_view(discretization),
      grid_transport_view(cell_transport_views),
      q_moments(source_moments),
      groupset(in_groupset),
      xsections(in_xsections),
      num_moms(in_num_moms),
      G(in_groupset.groups.size()),
      max_num_cell_dofs(in_max_num_cell_dofs),
      a_and_b_initialized(false),
      spls_index(0),
      angle_set_index(0)
  {}

  // ############################################################ Actual chunk
  void Sweep(chi_mesh::sweep_management::AngleSet* angle_set) override
  {
    if (!a_and_b_initialized)
    {
      Amat.resize(max_num_cell_dofs, std::vector<double>(max_num_cell_dofs));
      Atemp.resize(max_num_cell_dofs, std::vector<double>(max_num_cell_dofs));
      b.resize(G, std::vector<double>(max_num_cell_dofs, 0.0));
      source.resize(max_num_cell_dofs, 0.0);
      a_and_b_initialized = true;
    }

    const auto spds = angle_set->GetSPDS();
    const auto fluds = angle_set->fluds;

    const GsSubSet& subset = groupset.grp_subsets[angle_set->ref_subset];
    const int gs_ss_size  = groupset.grp_subset_sizes[angle_set->ref_subset];
    const int gs_ss_begin = subset.first;
    const int gs_gi = groupset.groups[gs_ss_begin].id; // Groupset subset first group number

    int deploc_face_counter = -1;
    int preloc_face_counter = -1;

    auto const& d2m_op = groupset.quadrature->GetDiscreteToMomentOperator();
    auto const& m2d_op = groupset.quadrature->GetMomentToDiscreteOperator();

    // ========================================================== Loop over each cell
    size_t num_loc_cells = spds->spls.item_id.size();
    for (spls_index = 0; spls_index < num_loc_cells; ++spls_index)
    {
      const int cell_local_id = spds->spls.item_id[spls_index];
      const auto& cell = grid_view->local_cells[cell_local_id];
      const auto& fe_intgrl_values = grid_fe_view.GetUnitIntegrals(cell);
      const int num_faces = cell.faces.size();
      const int num_dofs = fe_intgrl_values.NumNodes();
      const auto & transport_view = grid_transport_view[cell.local_id];
      const auto & sigma_tg = xsections[transport_view.xs_id]->sigma_tg;
      std::vector<bool> face_incident_flags(num_faces, false);

      // =================================================== Get Cell matrices
      const std::vector<std::vector<chi_mesh::Vector3>>& L =
        fe_intgrl_values.GetIntV_shapeI_gradshapeJ();

      const std::vector<std::vector<double>>& M =
        fe_intgrl_values.GetIntV_shapeI_shapeJ();

      const std::vector<std::vector<std::vector<double>>>& N =
        fe_intgrl_values.GetIntS_shapeI_shapeJ();

      // =================================================== Loop over angles in set
      const int ni_deploc_face_counter = deploc_face_counter;
      const int ni_preloc_face_counter = preloc_face_counter;
      for (angle_set_index = 0; angle_set_index < angle_set->angles.size(); ++angle_set_index)
      {
        deploc_face_counter = ni_deploc_face_counter;
        preloc_face_counter = ni_preloc_face_counter;
        const int angle_num = angle_set->angles[angle_set_index];
        const chi_mesh::Vector3 omega = groupset.quadrature->omegas[angle_num];

        // ============================================ Gradient matrix
        for (int i = 0; i < num_dofs; ++i)
          for (int j = 0; j < num_dofs; ++j)
            Amat[i][j] = omega.Dot(L[i][j]);

        for (int gsg = 0; gsg < gs_ss_size; ++gsg)
          b[gsg].assign(num_dofs, 0.0);

        // ============================================ Surface integrals
        int in_face_counter = -1;
        for (int f = 0; f < num_faces; ++f)
        {
          const auto& face = cell.faces[f];
          const double mu = omega.Dot(face.normal);

          if (mu < 0.0) // Upwind
          {
            face_incident_flags[f] = true;
            const bool local = transport_view.face_local[f];
            const bool boundary = not face.has_neighbor;
            const int num_face_indices = face.vertex_ids.size();
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
                  const double mu_Nij = -mu*N[f][i][j];
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
                  const double mu_Nij = -mu*N[f][i][j];
                  Amat[i][j] += mu_Nij;
                  for (int gsg = 0; gsg < gs_ss_size; ++gsg)
                    b[gsg][i] += psi[gsg]*mu_Nij;
                }
              }
            }
            else
            {
              // This counter update-logic is for mapping an incident boundary
              // condition. Because it is cheap, the cell faces was mapped to a
              // corresponding boundary during initialization and is
              // independent of angle. Accessing things like reflective boundary
              // angular fluxes (and complex boundary conditions), requires the
              // more general bndry_face_counter.
              const int bndry_index = face.neighbor_id;
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
                                                    suppress_surface_src);
                  const double mu_Nij = -mu*N[f][i][j];
                  Amat[i][j] += mu_Nij;
                  for (int gsg = 0; gsg < gs_ss_size; ++gsg)
                    b[gsg][i] += psi[gsg]*mu_Nij;
                }
              }
            }
          } // if upwind
        } // for f

        // ========================================== Looping over groups
        for (int gsg = 0; gsg < gs_ss_size; ++gsg)
        {
          const int g = gs_gi+gsg;

          // ============================= Contribute source moments
          for (int i = 0; i < num_dofs; ++i)
          {
            double temp_src = 0.0;
            for (int m = 0; m < num_moms; ++m)
            {
              const int ir = transport_view.MapDOF(i, m, g);
              temp_src += m2d_op[m][angle_num]*(*q_moments)[ir];
            }
            source[i] = temp_src;
          }

          // ============================= Mass Matrix and Source
          const double sigma_tgr = sigma_tg[g];
          for (int i = 0; i < num_dofs; ++i)
          {
            double temp = 0.0;
            for (int j = 0; j < num_dofs; ++j)
            {
              const double Mij = M[i][j];
              Atemp[i][j] = Amat[i][j] + Mij*sigma_tgr;
              temp += Mij*source[j];
            }
            b[gsg][i] += temp;
          }

          // ============================= Solve system
          chi_math::GaussElimination(Atemp, b[gsg], fe_intgrl_values.NumNodes());
        }

        // ============================= Accumulate flux
        for (int m = 0; m < num_moms; ++m)
        {
          const double wn_d2m = d2m_op[m][angle_num];
          for (int i = 0; i < num_dofs; ++i)
          {
            const int ir = transport_view.MapDOF(i, m, gs_gi);
            for (int gsg = 0; gsg < gs_ss_size; ++gsg)
              (*x)[ir + gsg] += wn_d2m*b[gsg][i];
          }
        }

        for (auto callback : moment_callbacks)
          callback(this, angle_set);

        // ============================= Save angular fluxes if needed
        if (groupset.psi_to_be_saved)
        {
          auto& psi = groupset.psi_new_local;
          auto& psi_uk_man = groupset.psi_uk_man;
          for (int i = 0; i < num_dofs; ++i)
          {
            int ir = grid_fe_view.MapDOFLocal(cell,i,psi_uk_man,angle_num,0);
            for (int gsg = 0; gsg < gs_ss_size; ++gsg)
              psi[ir + gsg] = b[gsg][i];
          }//for i
        }//if save psi

        int out_face_counter = -1;
        for (int f = 0; f < num_faces; ++f)
        {
          if (face_incident_flags[f]) continue;

          // ============================= Set flags and counters
          out_face_counter++;
          const auto& face = cell.faces[f];
          const bool local = transport_view.face_local[f];
          const bool boundary = not face.has_neighbor;
          const int num_face_indices = face.vertex_ids.size();

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
            const int bndry_index = face.neighbor_id;
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
      } // for n
    } // for cell
  }//Sweep
};

#endif
