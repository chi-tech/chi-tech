#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiPhysics/chi_physics.h"
#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Initializes parallel arrays.*/
void LinearBoltzmann::Solver::InitializeParrays()
{
  //================================================== Initialize unknown structure
  for (int m=0; m<num_moments; m++)
  {
    flux_moments_uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, groups.size());
    flux_moments_uk_man.unknowns.back().text_name = "m"+std::to_string(m);
  }

  //================================================== Compute local # of dof
  auto& per_node = ChiMath::UNITARY_UNKNOWN_MANAGER;
  local_node_count = discretization->GetNumLocalDOFs(per_node);
  glob_node_count = discretization->GetNumGlobalDOFs(per_node);

  //================================================== Compute num of unknowns
  size_t num_grps = groups.size();
  size_t local_unknown_count = local_node_count * num_grps * num_moments;

  chi_log.Log(LOG_ALLVERBOSE_1) << "LBS Number of phi unknowns: "
                                << local_unknown_count;

  //================================================== Size local vectors
  q_moments_local.assign(local_unknown_count,0.0);
  phi_old_local.assign(local_unknown_count,0.0);
  phi_new_local.assign(local_unknown_count,0.0);

  //============================================= Setup groupset psi vectors
  for (auto& groupset : groupsets)
  {
    psi_new_local.emplace_back();
    if (options.save_angular_flux)
    {
      size_t num_ang_unknowns =
          discretization->GetNumLocalDOFs(groupset.psi_uk_man);
      psi_new_local.back().assign(num_ang_unknowns, 0.0);
    }
  }

  //============================================= Setup precursor vector
  if (options.use_precursors)
  {
    size_t num_precursor_dofs =
        grid->local_cells.size() * max_precursors_per_material;
    precursor_new_local.assign(num_precursor_dofs, 0.0);
  }

  //================================================== Read Restart data
  if (options.read_restart_data)
    ReadRestartData(options.read_restart_folder_name,
                    options.read_restart_file_base);
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

  chi_mesh::Vector3 ihat(1.0, 0.0, 0.0);
  chi_mesh::Vector3 jhat(0.0, 1.0, 0.0);
  chi_mesh::Vector3 khat(0.0, 0.0, 1.0);

  auto pwl =
      std::dynamic_pointer_cast<SpatialDiscretization_FE>(discretization);

  cell_transport_views.clear();
  cell_transport_views.reserve(grid->local_cells.size());
  for (auto& cell : grid->local_cells)
  {
    size_t num_nodes  = discretization->GetCellNumNodes(cell);
    int    mat_id     = cell.material_id;
    int    xs_mapping = matid_to_xs_map[mat_id];

    //compute cell volumes
    double cell_volume = 0.0;
    auto& fe_values = pwl->GetUnitIntegrals(cell);
    for (int i = 0; i < fe_values.NumNodes(); ++i)
      cell_volume += fe_values.IntV_shapeI(i);

    size_t cell_phi_address = block_MG_counter;

    std::vector<bool> face_local_flags;
    face_local_flags.resize(cell.faces.size(), true);
    bool cell_on_boundary = false;
    int f=0;
    for (auto& face : cell.faces)
    {
      if (not face.has_neighbor)
      {
        chi_mesh::Vector3& n = face.normal;

        int boundary_id = -1;
        if      (n.Dot(ihat)>0.999)  boundary_id = 0;
        else if (n.Dot(ihat)<-0.999) boundary_id = 1;
        else if (n.Dot(jhat)> 0.999) boundary_id = 2;
        else if (n.Dot(jhat)<-0.999) boundary_id = 3;
        else if (n.Dot(khat)> 0.999) boundary_id = 4;
        else if (n.Dot(khat)<-0.999) boundary_id = 5;

        if (boundary_id >= 0) face.neighbor_id = boundary_id;
        cell_on_boundary = true;
      }//if bndry

      if (not face.IsNeighborLocal(*grid))
        face_local_flags[f] = false;
      ++f;
    }//for f

    if (num_nodes > max_cell_dof_count)
      max_cell_dof_count = num_nodes;

    cell_transport_views.emplace_back(cell_phi_address,
                                      num_nodes,
                                      num_grps,
                                      num_moments,
                                      xs_mapping,
                                      cell_volume,
                                      face_local_flags,
                                      cell_on_boundary);
    block_MG_counter += num_nodes * num_grps * num_moments;
  }//for local cell

  //================================================== Populate grid nodal mappings
  // This is used in the Flux Data Structures (FLUDS)
  grid_nodal_mappings.clear();
  grid_nodal_mappings.reserve(grid->local_cells.size());
  for (auto& cell : grid->local_cells)
  {
    chi_mesh::sweep_management::CellFaceNodalMapping cell_nodal_mapping;
    cell_nodal_mapping.reserve(cell.faces.size());

    for (auto& face : cell.faces)
    {
      std::vector<short> face_nodal_mapping;
      int ass_face = -1;

      if (face.has_neighbor and face.IsNeighborLocal(*grid))
      {
        grid->FindAssociatedVertices(face,face_nodal_mapping);
        ass_face = face.GetNeighborAssociatedFace(*grid);
      }

      cell_nodal_mapping.emplace_back(ass_face,face_nodal_mapping);
    }//for f

    grid_nodal_mappings.push_back(cell_nodal_mapping);
  }//for local cell

  //================================================== Initialize Field Functions
  if (field_functions.empty())
  {
    for (int g=0; g<groups.size(); g++)
    {
      for (int m=0; m<num_moments; m++)
      {
        std::string text_name = std::string("Flux_g") +
                                std::to_string(g) +
                                std::string("_m") + std::to_string(m);

        auto group_ff = std::make_shared<chi_physics::FieldFunction>(
          text_name,              //Field name
          discretization,         //Spatial discretization
          &phi_old_local,         //Data vector
          flux_moments_uk_man,    //Unknown manager
          m,                      //Reference unknown
          g);                     //Reference component

        chi_physics_handler.fieldfunc_stack.push_back(group_ff);
        field_functions.push_back(group_ff);
      }//for m
    }//for g
  }//if empty

}
