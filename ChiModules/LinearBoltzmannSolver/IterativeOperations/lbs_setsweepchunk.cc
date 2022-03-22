#include "../lbs_linear_boltzmann_solver.h"
#include "../SweepChunks/lbs_sweepchunk_pwl.h"

typedef chi_mesh::sweep_management::SweepChunk SweepChunk;

//###################################################################
/**Sets up the sweek chunk for the given discretization method.*/
std::shared_ptr<SweepChunk> lbs::SteadySolver::
  SetSweepChunk(LBSGroupset& groupset)
{
  auto pwl_sdm =
    std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(discretization);

  //================================================== Setting up required
  //                                                   sweep chunks
  auto sweep_chunk = std::make_shared<SweepChunkPWL>(
    grid,                                    //Spatial grid of cells
        *pwl_sdm,                                //Spatial discretization
        cell_transport_views,                    //Cell transport views
        phi_new_local,                           //Destination phi
        psi_new_local[groupset.id],              //Destination psi
        q_moments_local,                         //Source moments
        groupset,                                //Reference groupset
        matid_to_xs_map,                             //Material cross-sections
        num_moments,
    max_cell_dof_count);

  return sweep_chunk;
}
