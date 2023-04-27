#include "pi_keigen_scdsa.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include "ChiMath/VectorGhostCommunicator/vector_ghost_communicator.h"

namespace lbs
{

// ##################################################################
/**Copies only the scalar moments from an lbs primary flux moments
 * vector.*/
std::vector<double>
XXPowerIterationKEigenSCDSA::CopyOnlyPhi0(const LBSGroupset& groupset,
                                          const std::vector<double>& phi_in)
{
  typedef const int64_t cint64;

  const auto& lbs_sdm = lbs_solver_.SpatialDiscretization();
  const auto& diff_sdm = diffusion_solver_->SpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();
  const auto& phi_uk_man = lbs_solver_.UnknownManager();

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  const size_t diff_num_local_dofs =
    requires_ghosts_ ? diff_sdm.GetNumLocalAndGhostDOFs(diff_uk_man)
                     : diff_sdm.GetNumLocalDOFs(diff_uk_man);

  std::vector<double> phi_data =
    (not continuous_sdm_ptr_) ? phi_in : MakeContinuousVersion(phi_in);

  VecDbl output_phi_local(diff_num_local_dofs, 0.0);

  for (const auto& cell : lbs_solver_.Grid().local_cells)
  {
    const auto& cell_mapping = lbs_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; i++)
    {
      cint64 diff_phi_map = diff_sdm.MapDOFLocal(cell, i, diff_uk_man, 0, 0);
      cint64 lbs_phi_map = lbs_sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      double* output_mapped = &output_phi_local[diff_phi_map];
      const double* phi_in_mapped = &phi_data[lbs_phi_map];

      for (size_t g = 0; g < gss; g++)
      {
        output_mapped[g] = phi_in_mapped[g];
      } // for g
    }   // for node
  }     // for cell

  return output_phi_local;
}

// ##################################################################
/**Copies back only the scalar moments to a lbs primary flux vector.*/
void XXPowerIterationKEigenSCDSA::ProjectBackPhi0(
  const LBSGroupset& groupset,
  const std::vector<double>& input,
  std::vector<double>& output)
{
  typedef const int64_t cint64;

  const auto& lbs_sdm = lbs_solver_.SpatialDiscretization();
  const auto& diff_sdm = diffusion_solver_->SpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();
  const auto& phi_uk_man = lbs_solver_.UnknownManager();

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  const size_t diff_num_local_dofs =
    requires_ghosts_ ? diff_sdm.GetNumLocalAndGhostDOFs(diff_uk_man)
                     : diff_sdm.GetNumLocalDOFs(diff_uk_man);

  ChiLogicalErrorIf(input.size() != diff_num_local_dofs,
                    "Vector size mismatch");

  for (const auto& cell : lbs_solver_.Grid().local_cells)
  {
    const auto& cell_mapping = lbs_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; i++)
    {
      cint64 diff_phi_map = diff_sdm.MapDOFLocal(cell, i, diff_uk_man, 0, 0);
      cint64 lbs_phi_map = lbs_sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      const double* input_mapped = &input[diff_phi_map];
      double* output_mapped = &output[lbs_phi_map];

      for (int g = 0; g < gss; g++)
        output_mapped[g] = input_mapped[g];
    } // for dof
  }   // for cell
}

// ##################################################################
/**Creates a ghost communicator and all associated information.*/
XXPowerIterationKEigenSCDSA::GhostInfo
XXPowerIterationKEigenSCDSA::MakePWLDVecGhostCommInfo()
{
  chi::log.Log() << "Making PWLD ghost communicator";
  const auto& lbs_uk_man = lbs_solver_.UnknownManager();
  const auto& lbs_sdm = lbs_solver_.SpatialDiscretization();

  const size_t num_local_dofs = lbs_sdm.GetNumLocalDOFs(lbs_uk_man);
  const size_t num_globl_dofs = lbs_sdm.GetNumGlobalDOFs(lbs_uk_man);

  chi::log.Log() << "Number of global dofs" << num_globl_dofs;

  const size_t num_unknowns = lbs_uk_man.unknowns_.size();

  // Build a list of global ids
  std::set<int64_t> global_dof_ids_set;

  const auto& grid = lbs_solver_.Grid();
  const auto ghost_cell_ids = grid.cells.GetGhostGlobalIDs();
  for (const auto global_id : ghost_cell_ids)
  {
    const auto& cell = grid.cells[global_id];
    const auto& cell_mapping = lbs_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t u = 0; u < num_unknowns; ++u)
      {
        const size_t num_comps = lbs_uk_man.unknowns_[u].num_components_;
        for (size_t c = 0; c < num_comps; ++c)
        {
          const int64_t dof_map = lbs_sdm.MapDOF(cell, i, lbs_uk_man, u, c);
          global_dof_ids_set.insert(dof_map);
        } // for component
      }   // for unknown
    }     // for node i
  }       // for ghost cell

  // Convert the list to a vector
  std::vector<int64_t> global_indices(global_dof_ids_set.begin(),
                                      global_dof_ids_set.end());

  // Create the vector ghost communicator
  auto vgc = std::make_shared<chi_math::VectorGhostCommunicator>(
    num_local_dofs, num_globl_dofs, global_indices, MPI_COMM_WORLD);

  // Create the map
  std::map<int64_t, int64_t> ghost_global_id_2_local_map;
  {
    int64_t k = 0;
    for (const auto ghost_id : global_indices)
    {
      ghost_global_id_2_local_map[ghost_id] =
        static_cast<int64_t>(num_local_dofs + k++);
    }
  }

  chi::log.Log() << "Done making PWLD ghost communicator";
  return {vgc, ghost_global_id_2_local_map};
}

// ##################################################################
/**Makes a continuous version of a pwld lbs primary flux vector that
 * is still technically in pwld format but the local nodes have been
 * averaged at the discontinuities.*/
std::vector<double> XXPowerIterationKEigenSCDSA::MakeContinuousVersion(
  const std::vector<double>& input)
{
  ChiInvalidArgumentIf(
    not continuous_sdm_ptr_,
    "Requires a continuous spatial discretization to be defined.");

  ChiLogicalErrorIf(not pwld_ghost_info_.vector_ghost_communicator,
                    "No vector ghost communicator defined for pwld_ghost_info");

  const auto& vgc = pwld_ghost_info_.vector_ghost_communicator;
  const auto& dfem_dof_global2local_map =
    pwld_ghost_info_.ghost_global_id_2_local_map;

  auto input_with_ghosts = vgc->MakeGhostedVector(input);
  vgc->CommunicateGhostEntries(input_with_ghosts);

  typedef const int64_t cint64_t;

  const auto& pwlc_sdm = *continuous_sdm_ptr_;
  const auto& pwld_sdm = lbs_solver_.SpatialDiscretization();
  const auto& uk_man = lbs_solver_.UnknownManager();
  const auto& grid = lbs_solver_.Grid();

  const size_t num_unknowns = uk_man.unknowns_.size();

  const size_t num_cfem_local_dofs = pwlc_sdm.GetNumLocalAndGhostDOFs(uk_man);

  std::vector<double> cont_input(num_cfem_local_dofs, 0.0);
  std::vector<double> cont_input_ctr(num_cfem_local_dofs, 0.0);

  std::map<int64_t, int64_t> cfem_dof_global2local_map;

  //================================================== Local cells first
  std::set<uint64_t> partition_bndry_vertex_id_set;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = pwld_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t u = 0; u < num_unknowns; ++u)
      {
        const size_t num_components = uk_man.unknowns_[u].num_components_;
        for (size_t c = 0; c < num_components; ++c)
        {
          cint64_t dof_dfem_map = pwld_sdm.MapDOFLocal(cell, i, uk_man, u, c);
          cint64_t dof_cfem_map = pwlc_sdm.MapDOFLocal(cell, i, uk_man, u, c);
          cint64_t dof_cfem_map_globl = pwlc_sdm.MapDOF(cell, i, uk_man, u, c);

          cfem_dof_global2local_map[dof_cfem_map_globl] = dof_cfem_map;

          const double phi_value = input[dof_dfem_map];

          cont_input[dof_cfem_map] += phi_value;
          cont_input_ctr[dof_cfem_map] += 1.0;
        } // for component c
      }   // for unknown u
    }     // for node i

    for (const auto& face : cell.faces_)
      if (face.has_neighbor_)
        if (not grid.IsCellLocal(face.neighbor_id_))
          for (const uint64_t vid : face.vertex_ids_)
            partition_bndry_vertex_id_set.insert(vid);
  } // for local cell

  //================================================== Ghost cells
  const auto ghost_cell_ids = grid.cells.GetGhostGlobalIDs();
  const auto& vid_set = partition_bndry_vertex_id_set;
  for (const auto global_id : ghost_cell_ids)
  {
    const auto& cell = grid.cells[global_id];
    const auto& cell_mapping = pwld_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      if (vid_set.find(cell.vertex_ids_[i]) == vid_set.end()) continue;

      for (size_t u = 0; u < num_unknowns; ++u)
      {
        const size_t num_components = uk_man.unknowns_[u].num_components_;
        for (size_t c = 0; c < num_components; ++c)
        {
          cint64_t dof_dfem_map_globl = pwld_sdm.MapDOF(cell, i, uk_man, u, c);
          cint64_t dof_cfem_map_globl = pwlc_sdm.MapDOF(cell, i, uk_man, u, c);
          if (cfem_dof_global2local_map.count(dof_cfem_map_globl) > 0)
          {
            cint64_t dof_dfem_map =
              dfem_dof_global2local_map.at(dof_dfem_map_globl);
            cint64_t dof_cfem_map =
              cfem_dof_global2local_map[dof_cfem_map_globl];

            const double phi_value = input_with_ghosts[dof_dfem_map];

            cont_input[dof_cfem_map] += phi_value;
            cont_input_ctr[dof_cfem_map] += 1.0;
          }
        } // for component
      }   // for unknown
    }     // for node i
  }       // for ghost cell

  //================================================== Compute nodal averages
  {
    const size_t num_vals = cont_input.size();
    for (size_t k = 0; k < num_vals; ++k)
      cont_input[k] /= cont_input_ctr[k];
  }

  //================================================== Project back to dfem
  auto output = input;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = pwld_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t u = 0; u < num_unknowns; ++u)
      {
        const size_t num_components = uk_man.unknowns_[u].num_components_;
        for (size_t c = 0; c < num_components; ++c)
        {
          cint64_t dof_dfem_map = pwld_sdm.MapDOFLocal(cell, i, uk_man, u, c);
          cint64_t dof_cfem_map = pwlc_sdm.MapDOFLocal(cell, i, uk_man, u, c);

          const double phi_value = cont_input[dof_cfem_map];

          output[dof_dfem_map] = phi_value;
        } // for component c
      }   // for unknown u
    }     // for node i
  }       // for local cell

  return output;
}

} // namespace lbs