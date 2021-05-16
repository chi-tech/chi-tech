#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiPhysics/chi_physics.h"
#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Initializes p_arrays.\n
The question arises of what datatype can store the total amount of
unknowns. For now we will say we want to be
designing for 100 billion cells with
an assumed shape of a truncated octahedron which has 24 vertices.
We will also assume that we will be able to do 2000 energy groups
and finally we will assume we will do scattering orders up to 16
which requires 289 moments.
   DOFS per truncated octahedron = 24\n
   Energy groups                 = 2000\n
   Moments                       = 289\n
   # of cells                    =   100,000,000,000\n
   Total DOFS                    = 2,400,000,000,000\n
   Unknowns per cell             =        13,872,000\n
   Total Unknowns                = A crap ton\n
\n
It is easy to see here that this is a hell of a lot so how about we think about
something more modest. Like 200 energy groups scattering order 5 (36 moments)
and 2 billion cells.\n
   Energy groups                 = 200\n
   Moments                       = 36\n
   # of cells                    =     2,000,000,000\n
   Total DOFS                    =    48,000,000,000\n
   Unknowns per cell             =             7,200\n
   Total Unknowns                = 1.44xe13\n
\n
A long int only supports up to 4.29e9. This obviously requires
unsigned long long int which can hold up to 2x2e63.\n
\n
Another interesting aspect is what it will take to get to exascale. For a
discrete ordinates code this will undoubtly be evident in the amount of angular flux
unknowns. 1 billion cells, 24 vertices, 200 groups, 48 azimuthal angles per
octant, 8 polar angles per octant (3072) angles. 1.47456e16. Just a factor 68
away from exascale.
   */
void LinearBoltzmann::Solver::InitializeParrays()
{
  auto pwl_discretization =
    std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(discretization);

  //================================================== Initialize unknown structure
  for (int m=0; m<num_moments; m++)
  {
    flux_moments_uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, groups.size());
    auto& moment = flux_moments_uk_man.unknowns.back().text_name = "m"+std::to_string(m);
  }

  //================================================== Compute local # of dof
  auto GxM = flux_moments_uk_man.GetTotalUnknownStructureSize();
  local_node_count = pwl_discretization->GetNumLocalDOFs(flux_moments_uk_man) / GxM;
  globl_node_count = pwl_discretization->GetNumGlobalDOFs(flux_moments_uk_man) / GxM;

  //================================================== Compute num of unknowns
  size_t num_grps = groups.size();
  size_t M = num_moments;
  size_t local_dof_count = local_node_count * num_grps * M;

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
                                      num_grps, M);
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
