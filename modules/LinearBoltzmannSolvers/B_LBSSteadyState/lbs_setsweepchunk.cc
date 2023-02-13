#include "B_LBSSteadyState/lbs_linear_boltzmann_solver.h"
#include "B_LBSSteadyState/SweepChunks/lbs_sweepchunk_pwl.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"


typedef chi_mesh::sweep_management::SweepChunk SweepChunk;

//###################################################################
/**Sets up the sweek chunk for the given discretization_ method.*/
std::shared_ptr<SweepChunk> lbs::SteadyStateSolver::
  SetSweepChunk(LBSGroupset& groupset)
{
  auto pwl_sdm =
    std::dynamic_pointer_cast<chi_math::SpatialDiscretization_PWLD>(discretization_);

  //================================================== Setting up required
  //                                                   sweep chunks
  auto sweep_chunk = std::make_shared<SweepChunkPWL>(
    grid_ptr_,                                    //Spatial grid_ptr_ of cells
        *pwl_sdm,                                //Spatial discretization_
        cell_transport_views_,                    //Cell transport views
        phi_new_local_,                           //Destination phi
        psi_new_local_[groupset.id],              //Destination psi
        q_moments_local_,                         //Source moments
        groupset,                                //Reference groupset
        matid_to_xs_map_,                             //Material cross-sections
        num_moments_,
    max_cell_dof_count_);

  return sweep_chunk;
}
