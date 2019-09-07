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
  if (chi_physics_handler.material_stack.size() == 0)
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
    for (int p=0; p<cur_mat->properties.size(); p++)
    {
      if (cur_mat->properties[p]->type_index == TRANSPORT_XSECTIONS)
      {
        chi_physics::TransportCrossSections* transp_xs =
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



  //=================================== Initialize tallies
  size_t tally_size = num_grps*grid->local_cell_glob_indices.size();
  phi_tally_contrib.resize(tally_size,0.0);
  phi_tally.resize(tally_size,0.0);
  phi_tally_sqr.resize(tally_size,0.0);

  phi_global.resize(tally_size,0.0);
  phi_global_tally_sqr.resize(tally_size,0.0);

  phi_local_relsigma.resize(tally_size,0.0);

  //=================================== Initialize discretization
  fv_discretization = new SpatialDiscretization_FV;

  fv_discretization->
  AddViewOfLocalContinuum(grid,
                          grid->local_cell_glob_indices.size(),
                          grid->local_cell_glob_indices.data());

  //=================================== Initialize Sources
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



  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}