#include"solver_montecarlon.h"
#include <ChiPhysics/chi_physics.h>
#include <ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h>
#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>

#include <ChiMath/Statistics/cdfsampler.h>

extern ChiPhysics chi_physics_handler;

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog chi_log;
extern ChiMPI chi_mpi;

typedef unsigned long long TULL;

#include <ChiTimer/chi_timer.h>

/**Initialize the solver*/
bool chi_montecarlon::Solver::Initialize()
{
  chi_log.Log(LOG_0) << "Initializing MonteCarlo solver.";

  //=================================== Initialize Materials
  if (chi_physics_handler.material_stack.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_montecarlon::Solver::Initialize : No materials found.";
    exit(EXIT_FAILURE);
  }
  int num_mat = chi_physics_handler.material_stack.size();
  matid_xs_map.resize(num_mat,-1);
  for (int m=0; m<num_mat; m++)
  {
    chi_physics::Material* cur_mat = chi_physics_handler.material_stack[m];

    //======================= Only first xs will be used
    size_t num_props = cur_mat->properties.size();
    for (size_t p=0; p<num_props; p++)
    {
      if (cur_mat->properties[p]->type_index == TRANSPORT_XSECTIONS)
      {
        auto transp_xs =
          (chi_physics::TransportCrossSections*)cur_mat->properties[p];

        transp_xs->ComputeDiffusionParameters();
        transp_xs->ComputeDiscreteScattering(scattering_order);

        if (transp_xs->G > num_grps)
          num_grps = transp_xs->G;

        matid_xs_map[m] = p;
        break;
      }
    }//for prop
  }//for mat
  MPI_Barrier(MPI_COMM_WORLD);

  chi_mesh::Region*  aregion = this->regions.back();
  this->grid                 = aregion->volume_mesh_continua.back();



  //=================================== Process rendesvous intervals
  int num_batches = std::ceil((double)num_particles/tally_rendezvous_intvl);
  batch_sizes.resize(num_batches,tally_rendezvous_intvl);
  batch_sizes[num_batches-1] = num_particles -
                              (num_batches-1)*tally_rendezvous_intvl;

  chi_log.Log(LOG_0) << "Number of MPI merges to be performed: " << num_batches;

  //=================================== Procces num_part per interval
  for (int b=0; b<num_batches; b++)
  {
    TULL loc_num_part = std::ceil((double)batch_sizes[b]/chi_mpi.process_count);

    batch_sizes_per_loc.push_back(loc_num_part);
    if (chi_mpi.location_id == (chi_mpi.process_count-1))
    {
      batch_sizes_per_loc[b] = batch_sizes[b] -
                               chi_mpi.location_id*loc_num_part;
    }
  }
  chi_log.Log(LOG_0) << "Batches seperated among processes.";
  MPI_Barrier(MPI_COMM_WORLD);

  //=================================== Initialize tallies
  size_t num_local_cells = grid->local_cell_glob_indices.size();
  size_t tally_size = num_grps*num_local_cells;
  phi_tally_contrib.resize(tally_size,0.0);
  phi_tally.resize(tally_size,0.0);
  phi_tally_sqr.resize(tally_size,0.0);

  phi_global.resize(tally_size,0.0);
  phi_global_tally_sqr.resize(tally_size,0.0);

  phi_local_relsigma.resize(tally_size,0.0);

  //=================================== Initialize default discretization
  chi_log.Log(LOG_0) << "Adding finite volume views.";
  fv_discretization = new SpatialDiscretization_FV;

  fv_discretization->
    AddViewOfLocalContinuum(grid,
                            num_local_cells,
                            grid->local_cell_glob_indices.data());
  MPI_Barrier(MPI_COMM_WORLD);
  if (make_pwld)
  {
    chi_log.Log(LOG_0) << "Adding PWL finite element views.";
    //=================================== Initialize pwl discretization
    pwl_discretization = new SpatialDiscretization_PWL;

    pwl_discretization->
      AddViewOfLocalContinuum(grid,
                              num_local_cells,
                              grid->local_cell_glob_indices.data());

    //=================================== Generate moment wise addresses
    num_moms = 1;
    local_cell_pwl_dof_array_address.resize(num_local_cells,0);
    int block_MG_counter = 0;
    for (size_t lc=0; lc<num_local_cells; lc++)
    {
      local_cell_pwl_dof_array_address[lc] = block_MG_counter;
      int cell_g_index = grid->local_cell_glob_indices[lc];
      auto cell_pwl_view = pwl_discretization->MapFeView(cell_g_index);

      block_MG_counter += cell_pwl_view->dofs*num_grps*num_moms;
    }

    //=================================== Initialize PWLD tallies
    tally_size = block_MG_counter-1;
    phi_pwl_tally_contrib.resize(tally_size,0.0);
    phi_pwl_tally.resize(tally_size,0.0);
    phi_pwl_tally_sqr.resize(tally_size,0.0);

    phi_pwl_global.resize(tally_size,0.0);
    phi_pwl_global_tally_sqr.resize(tally_size,0.0);

    phi_pwl_local_relsigma.resize(tally_size,0.0);
  }
  MPI_Barrier(MPI_COMM_WORLD);


  //=================================== Initialize Sources
  chi_log.Log(LOG_0) << "Initializing sources";
  for (int s=0; s<sources.size(); s++)
  {
    sources[s]->Initialize(grid,fv_discretization);
  }

  //=================================== Initialize field functions
  for (int g=0; g<num_grps; g++)
  {
    chi_physics::FieldFunction* group_ff =
      new chi_physics::FieldFunction;
    group_ff->text_name = std::string("Flux_g") +
                          std::to_string(g);

    group_ff->grid = grid;
    group_ff->spatial_discretization = fv_discretization;
    group_ff->id = chi_physics_handler.fieldfunc_stack.size();

    group_ff->type = FF_SDM_FV;
    group_ff->num_grps = num_grps;
    group_ff->num_moms = 1;
    group_ff->grp = g;
    group_ff->mom = 0;
    group_ff->field_vector_local = &phi_global;

    chi_physics_handler.fieldfunc_stack.push_back(group_ff);
    field_functions.push_back(group_ff);
  }

  if (make_pwld)
  {
    for (int g=0; g<num_grps; g++)
    {
      for (int m=0; m<num_moms; m++)
      {
        chi_physics::FieldFunction* group_ff =
          new chi_physics::FieldFunction;
        group_ff->text_name = std::string("Flux_g") +
                              std::to_string(g) +
                              std::string("_m") + std::to_string(m);
        group_ff->grid = grid;
        group_ff->spatial_discretization = pwl_discretization;
        group_ff->id = chi_physics_handler.fieldfunc_stack.size();

        group_ff->type = FF_SDM_PWLD;
        group_ff->num_grps = num_grps;
        group_ff->num_moms = num_moms;
        group_ff->grp = g;
        group_ff->mom = m;
        group_ff->field_vector_local = &phi_pwl_global;
        group_ff->local_cell_dof_array_address = &local_cell_pwl_dof_array_address;

        chi_physics_handler.fieldfunc_stack.push_back(group_ff);
        field_functions.push_back(group_ff);
      }//for m
    }//for g
  }//if make_pwld




  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}