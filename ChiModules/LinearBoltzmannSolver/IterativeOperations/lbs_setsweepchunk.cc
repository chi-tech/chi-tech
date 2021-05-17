#include "../lbs_linear_boltzmann_solver.h"
#include "../SweepChunks/lbs_sweepchunk_pwl.h"

typedef chi_mesh::sweep_management::SweepChunk SweepChunk;

//###################################################################
/**Sets up the sweek chunk for the given discretization method.*/
std::shared_ptr<SweepChunk> LinearBoltzmann::Solver::SetSweepChunk(LBSGroupset& groupset)
{
  auto pwl_sdm = std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(discretization);

  if (not pwl_sdm)
    throw std::logic_error(std::string(__FUNCTION__) +
                           ": Trouble obtaining pwl spatial discretization.");

  //================================================== Setting up required
  //                                                   sweep chunks
  auto sweep_chunk = std::make_shared<LBSSweepChunkPWL>(
    grid,                                    //Spatial grid of cells
    *pwl_sdm,                                //Spatial discretization
    cell_transport_views,                    //Cell transport views
    &phi_new_local,                          //Destination phi
    &q_moments_local,                        //Source moments
    groupset,                                //Reference groupset
    material_xs,                             //Material cross-sections
    static_cast<int>(num_moments),           //Number of moments
    static_cast<int>(max_cell_node_count));  //Max num nodes per cell

  return sweep_chunk;

  //More will be added here as we have different SDMs
}
