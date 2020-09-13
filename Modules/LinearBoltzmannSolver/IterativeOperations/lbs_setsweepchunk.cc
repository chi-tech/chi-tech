#include "../lbs_linear_boltzmann_solver.h"
#include "../SweepChunks/lbs_sweepchunk_pwl.h"

typedef chi_mesh::sweep_management::SweepChunk SweepChunk;

//###################################################################
/**Sets up the sweek chunk for the given discretization method.*/
SweepChunk* LinearBoltzmann::Solver::SetSweepChunk(int group_set_num)
{
  //================================================== Obtain groupset
  LBSGroupset* groupset = group_sets[group_set_num];

  //================================================== Setting up required
  //                                                   sweep chunks
  SweepChunk* sweep_chunk = new LBSSweepChunkPWL(
        grid,                                    //Spatial grid of cells
        (SpatialDiscretization_PWL*)discretization, //Spatial discretization
        &cell_transport_views,                   //Cell transport views
        *this,
        &phi_new_local,                          //Destination phi
        &q_moments_local,                        //Source moments
        groupset,                                //Reference groupset
        &material_xs,                            //Material cross-sections
        num_moments,max_cell_dof_count);

  return sweep_chunk;
}