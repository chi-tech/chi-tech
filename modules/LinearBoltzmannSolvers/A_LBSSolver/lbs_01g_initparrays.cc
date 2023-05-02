#include "lbs_solver.h"

#include "ChiPhysics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"
#include "ChiConsole/chi_console.h"


//###################################################################
/**Initializes parallel arrays.*/
void lbs::LBSSolver::InitializeParrays()
{
  chi::log.Log() << "Initializing parallel arrays."
                 << " G=" << num_groups_
                 << " M=" << num_moments_ << std::endl;

  //================================================== Initialize unknown structure
  flux_moments_uk_man_.unknowns_.clear();
  for (size_t m=0; m < num_moments_; m++)
  {
    flux_moments_uk_man_.AddUnknown(chi_math::UnknownType::VECTOR_N, groups_.size());
    flux_moments_uk_man_.unknowns_.back().text_name_ = "m" + std::to_string(m);
  }

  //================================================== Compute local # of dof
  auto per_node = chi_math::UnknownManager::GetUnitaryUnknownManager();
  local_node_count_ = discretization_->GetNumLocalDOFs(per_node);
  glob_node_count_ = discretization_->GetNumGlobalDOFs(per_node);

  //================================================== Compute num of unknowns
  size_t num_grps = groups_.size();
  size_t local_unknown_count = local_node_count_ * num_grps * num_moments_;

  chi::log.LogAllVerbose1() << "LBS Number of phi unknowns: "
                                << local_unknown_count;

  //================================================== Size local vectors
  q_moments_local_.assign(local_unknown_count, 0.0);
  phi_old_local_.assign(local_unknown_count, 0.0);
  phi_new_local_.assign(local_unknown_count, 0.0);

  //============================================= Setup groupset psi vectors
  psi_new_local_.clear();
  for (auto& groupset : groupsets_)
  {
    psi_new_local_.emplace_back();
    if (options_.save_angular_flux)
    {
      size_t num_ang_unknowns =
          discretization_->GetNumLocalDOFs(groupset.psi_uk_man_);
      psi_new_local_.back().assign(num_ang_unknowns, 0.0);
    }
  }

  //============================================= Setup precursor vector
  if (options_.use_precursors)
  {
    size_t num_precursor_dofs =
      grid_ptr_->local_cells.size() * max_precursors_per_material_;
    precursor_new_local_.assign(num_precursor_dofs, 0.0);
  }

  //================================================== Read Restart data
  if (options_.read_restart_data)
    ReadRestartData(options_.read_restart_folder_name,
                    options_.read_restart_file_base);
  MPI_Barrier(MPI_COMM_WORLD);

  //================================================== Initialize transport views
  // Transport views act as a data structure to store information
  // related to the transport simulation. The most prominent function
  // here is that it holds the means to know where a given cell's
  // transport quantities are located in the unknown vectors (i.e. phi)
  //
  // Also, for a given cell, within a given sweep chunk,
  // we need to solve a matrix which square size is the
  // amount of nodes on the cell. max_cell_dof_count is
  // initialized here.
  //
  size_t block_MG_counter = 0;       //Counts the strides of moment and group

  const chi_mesh::Vector3 ihat(1.0, 0.0, 0.0);
  const chi_mesh::Vector3 jhat(0.0, 1.0, 0.0);
  const chi_mesh::Vector3 khat(0.0, 0.0, 1.0);

  max_cell_dof_count_ = 0;
  cell_transport_views_.clear();
  cell_transport_views_.reserve(grid_ptr_->local_cells.size());
  for (auto& cell : grid_ptr_->local_cells)
  {
    size_t num_nodes  = discretization_->GetCellNumNodes(cell);
    int    mat_id     = cell.material_id_;

    //compute cell volumes
    double cell_volume = 0.0;
    const auto& IntV_shapeI = unit_cell_matrices_[cell.local_id_].Vi_vectors;
    for (size_t i = 0; i < num_nodes; ++i)
      cell_volume += IntV_shapeI[i];

    size_t cell_phi_address = block_MG_counter;

    std::vector<bool> face_local_flags;
    face_local_flags.resize(cell.faces_.size(), true);
    bool cell_on_boundary = false;
    int f=0;
    for (auto& face : cell.faces_)
    {
      if (not face.has_neighbor_)
      {
        chi_mesh::Vector3& n = face.normal_;

        int boundary_id = -1;
        if      (n.Dot(ihat)>0.999)  boundary_id = 0;
        else if (n.Dot(ihat)<-0.999) boundary_id = 1;
        else if (n.Dot(jhat)> 0.999) boundary_id = 2;
        else if (n.Dot(jhat)<-0.999) boundary_id = 3;
        else if (n.Dot(khat)> 0.999) boundary_id = 4;
        else if (n.Dot(khat)<-0.999) boundary_id = 5;

        if (boundary_id >= 0) face.neighbor_id_ = boundary_id;
        cell_on_boundary = true;
      }//if bndry

      if (not face.IsNeighborLocal(*grid_ptr_))
        face_local_flags[f] = false;
      ++f;
    }//for f

    if (num_nodes > max_cell_dof_count_)
      max_cell_dof_count_ = num_nodes;

    cell_transport_views_.emplace_back(cell_phi_address,
                                       num_nodes,
                                       num_grps,
                                       num_moments_,
                                       *matid_to_xs_map_[mat_id],
                                       cell_volume,
                                       face_local_flags,
                                       cell_on_boundary);
    block_MG_counter += num_nodes * num_grps * num_moments_;
  }//for local cell

  //================================================== Populate grid nodal mappings
  // This is used in the Flux Data Structures (FLUDS)
  grid_nodal_mappings_.clear();
  grid_nodal_mappings_.reserve(grid_ptr_->local_cells.size());
  for (auto& cell : grid_ptr_->local_cells)
  {
    chi_mesh::sweep_management::CellFaceNodalMapping cell_nodal_mapping;
    cell_nodal_mapping.reserve(cell.faces_.size());

    for (auto& face : cell.faces_)
    {
      std::vector<short> face_nodal_mapping;
      int ass_face = -1;

      if (face.has_neighbor_ and face.IsNeighborLocal(*grid_ptr_))
      {
        grid_ptr_->FindAssociatedVertices(face, face_nodal_mapping);
        ass_face = face.GetNeighborAssociatedFace(*grid_ptr_);
      }

      cell_nodal_mapping.emplace_back(ass_face,face_nodal_mapping);
    }//for f

    grid_nodal_mappings_.push_back(cell_nodal_mapping);
  }//for local cell

  //================================================== Get grid localized
  //                                                   communicator set
  grid_local_comm_set_ = grid_ptr_->MakeMPILocalCommunicatorSet();

  //================================================== Make face histogram
  grid_face_histogram_ = grid_ptr_->MakeGridFaceHistogram();

  //================================================== Initialize Field Functions
  if (field_functions_.empty())
  {
    for (size_t g = 0; g < groups_.size(); ++g)
    {
      for (size_t m=0; m < num_moments_; m++)
      {
        std::string solver_name;
        if (not TextName().empty()) solver_name = TextName();

        char buff[100];
        snprintf(buff, 99, "%s_Flux_g%03d_m%02d",
                 solver_name.c_str(),
                 static_cast<int>(g),
                 static_cast<int>(m));
        const std::string text_name = std::string(buff);

        using namespace chi_math;
        auto group_ff = std::make_shared<chi_physics::FieldFunctionGridBased>(
          text_name,                     //Field name
          discretization_,                //Spatial discretization
          Unknown(UnknownType::SCALAR)); //Unknown/Variable

        chi::field_function_stack.push_back(group_ff);
        field_functions_.push_back(group_ff);
      }//for m
    }//for g
  }//if empty

  MPI_Barrier(MPI_COMM_WORLD);
  chi::log.Log()
    << "Done with parallel arrays.                Process memory = "
    << std::setprecision(3)
    << chi_objects::ChiConsole::GetMemoryUsageInMB() << " MB" << std::endl;
}
