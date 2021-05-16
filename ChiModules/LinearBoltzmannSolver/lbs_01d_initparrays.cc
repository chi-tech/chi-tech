#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "chi_log.h"
#include "ChiPhysics/chi_physics.h"

//###################################################################
/**Initializes data arrays and other data.*/
void LinearBoltzmann::Solver::InitializeParrays()
{
  auto& chi_log = ChiLog::GetInstance();
  auto& physics_handler = ChiPhysics::GetInstance();

  auto pwl_discretization =
    std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(discretization);

  if (not pwl_discretization)
    throw std::logic_error(std::string(__FUNCTION__) +
                           ": Unknown trouble with spatial discretization.");

  //================================================== Initialize unknown structure
  for (int m=0; m<num_moments; m++)
  {
    flux_moments_uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, groups.size());
    flux_moments_uk_man.unknowns.back().text_name = "m"+std::to_string(m);
  }

  //================================================== Compute local # of dof
  auto& per_node = ChiMath::UNITARY_UNKNOWN_MANAGER;
  local_node_count = pwl_discretization->GetNumLocalDOFs(per_node);
  globl_node_count = pwl_discretization->GetNumGlobalDOFs(per_node);

  //================================================== Compute num of unknowns
  size_t num_grps = groups.size();
  size_t local_dof_count = local_node_count * num_grps * num_moments;

  chi_log.Log(LOG_ALLVERBOSE_1) << "LBS Number of phi unknowns: "
                                << local_dof_count;

  //================================================== Size local vectors
  q_moments_local.assign(local_dof_count, 0.0);
  phi_old_local.assign(local_dof_count, 0.0);
  phi_new_local.assign(local_dof_count, 0.0);

  //================================================== Read Restart data
  if (options.read_restart_data)
    ReadRestartData(options.read_restart_folder_name,
                    options.read_restart_file_base);
  MPI_Barrier(MPI_COMM_WORLD);

  //================================================== Initialize transport views
  // Transport-views act as a data structure to store information
  // related to the transport simulation. The most prominent function
  // here is that it holds the means to know where a given cell's
  // transport quantities are located in the unknown vectors (i.e. phi)
  //
  // Also, for a given cell, within a given sweep chunk,
  // we need to solve a matrix which square size is the
  // amount of nodes on the cell. max_cell_node_count is
  // initialized here.
  //
  size_t block_MG_counter = 0;       //Counts the strides of moment and group

  chi_mesh::Vector3 ihat(1.0, 0.0, 0.0);
  chi_mesh::Vector3 jhat(0.0, 1.0, 0.0);
  chi_mesh::Vector3 khat(0.0, 0.0, 1.0);

  cell_transport_views.clear();
  for (auto& cell : grid->local_cells)
  {
    const auto& fe_intgrl_values = pwl_discretization->GetUnitIntegrals(cell);
    size_t cell_num_nodes = fe_intgrl_values.NumNodes();

    size_t cell_phi_address = block_MG_counter;
    block_MG_counter += cell_num_nodes * num_grps * num_moments;

    //Init face upwind flags and adj_partition_id
    std::vector<bool> face_local_flags;
    face_local_flags.resize(cell.faces.size(), true);
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
      }//if bndry

      if (not face.IsNeighborLocal(*grid))
        face_local_flags[f] = false;

      ++f;
    }//for f

    max_cell_node_count = std::max(max_cell_node_count,cell_num_nodes);

    cell_transport_views.emplace_back(cell_phi_address,
                                      fe_intgrl_values.NumNodes(),
                                      matid_to_xs_map[cell.material_id],
                                      face_local_flags,
                                      num_grps, num_moments);
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

        physics_handler.fieldfunc_stack.push_back(group_ff);
        field_functions.push_back(group_ff);
      }//for m
    }//for g
  }//if empty

}
