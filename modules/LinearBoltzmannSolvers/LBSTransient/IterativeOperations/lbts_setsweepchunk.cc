#include "LBSTransient/lbts_transient_solver.h"
#include "LBSTransient/SweepChunks/lbts_sweepchunk_pwl.h"

//###################################################################
/**Sets up the sweek chunk for the given discretization method.*/
std::shared_ptr<SweepChunk> lbs::TransientSolver::
  SetTransientSweepChunk(LBSGroupset& groupset)
{
  auto pwl_sdm =
    std::dynamic_pointer_cast<chi_math::SpatialDiscretization_PWLD>(discretization);

  double theta;
  if (method == chi_math::SteppingMethod::BACKWARD_EULER)
    theta = 1.0;
  else
    theta = 0.5;

  //================================================== Setting up required
  //                                                   sweep chunks
  auto sweep_chunk = std::make_shared<SweepChunkPWLTransientTheta>(
    grid,                                    //Spatial grid of cells
    *pwl_sdm,                                //Spatial discretization
    cell_transport_views,                    //Cell transport views
    phi_new_local,                           //Destination phi
    psi_new_local[groupset.id],              //Destination psi

    psi_prev_local[groupset.id],
    theta,
    dt,

    q_moments_local,                         //Source moments
    groupset,                                //Reference groupset
    matid_to_xs_map,                             //Material cross-sections
    num_moments,
    max_cell_dof_count);

  return sweep_chunk;
}