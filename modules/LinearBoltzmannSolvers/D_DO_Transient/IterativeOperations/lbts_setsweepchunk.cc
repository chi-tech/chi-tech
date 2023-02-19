#include "D_DO_Transient/lbts_transient_solver.h"
#include "D_DO_Transient/SweepChunks/lbts_sweepchunk_pwl.h"

//###################################################################
/**Sets up the sweek chunk for the given discretization_ method.*/
std::shared_ptr<SweepChunk> lbs::DiscOrdTransientSolver::
  SetTransientSweepChunk(LBSGroupset& groupset)
{
  auto pwl_sdm =
    std::dynamic_pointer_cast<chi_math::SpatialDiscretization_PWLD>(discretization_);

  double theta;
  if (method == chi_math::SteppingMethod::BACKWARD_EULER)
    theta = 1.0;
  else
    theta = 0.5;

  //================================================== Setting up required
  //                                                   sweep chunks
  auto sweep_chunk = std::make_shared<SweepChunkPWLTransientTheta>(
    grid_ptr_,                                    //Spatial grid_ptr_ of cells
    *pwl_sdm,                                //Spatial discretization_
    cell_transport_views_,                    //Cell transport views
    phi_new_local_,                           //Destination phi
    psi_new_local_[groupset.id],              //Destination psi

    psi_prev_local_[groupset.id],
    theta,
    dt_,

    q_moments_local_,                         //Source moments
    groupset,                                //Reference groupset
    matid_to_xs_map_,                             //Material cross-sections
    num_moments_,
    max_cell_dof_count_);

  return sweep_chunk;
}