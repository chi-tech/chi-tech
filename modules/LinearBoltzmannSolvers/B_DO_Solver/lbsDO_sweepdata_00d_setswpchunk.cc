#include "lbs_discrete_ordinates_solver.h"

#include "SweepChunks/lbs_sweepchunk_pwl.h"

typedef chi_mesh::sweep_management::SweepChunk SweepChunk;

//###################################################################
/**Sets up the sweek chunk for the given discretization method.*/
std::shared_ptr<SweepChunk> lbs::DiscreteOrdinatesSolver::
  SetSweepChunk(LBSGroupset& groupset)
{
  //================================================== Setting up required
  //                                                   sweep chunks
  auto sweep_chunk = std::make_shared<SweepChunkPWL>(
    grid_ptr_,                                    //Spatial grid of cells
    *discretization_,                             //Spatial discretization
    unit_cell_matrices_,                          //Unit cell matrices
    cell_transport_views_,                        //Cell transport views
    phi_new_local_,                               //Destination phi
    psi_new_local_[groupset.id_],                 //Destination psi
    q_moments_local_,                             //Source moments
    groupset,                                     //Reference groupset
    matid_to_xs_map_,                             //Material cross-sections
    num_moments_,
    max_cell_dof_count_);

  return sweep_chunk;
}