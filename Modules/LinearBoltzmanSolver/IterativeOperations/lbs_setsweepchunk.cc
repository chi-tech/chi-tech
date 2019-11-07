#include "../lbs_linear_boltzman_solver.h"
#include "../SweepChunks/lbs_sweepchunk_pwl.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/chi_volumemesher.h>
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include <ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h>
#include <ChiMesh/VolumeMesher/Predefined2D/volmesher_predefined2d.h>

typedef chi_mesh::sweep_management::SweepChunk SweepChunk;

//###################################################################
/**Sets up the sweek chunk for the given discretization method.*/
SweepChunk* LinearBoltzman::Solver::SetSweepChunk(int group_set_num)
{
  //================================================== Obtain groupset
  LBSGroupset* groupset = group_sets[group_set_num];

  //================================================== Obtain the mesher
  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  //================================================== Setting up required
  //                                                   sweep chunks
  SweepChunk* sweep_chunk = nullptr;

  sweep_chunk =
      new LBSSweepChunkPWL(
        grid,                                    //Spatial grid of cells
        (SpatialDiscretization_PWL*)discretization, //Spatial discretization
        &cell_transport_views,                   //Cell transport views
        &phi_new_local,                          //Destination phi
        &q_moments_local,                        //Source moments
        groupset,                                //Reference groupset
        &material_xs,                            //Material cross-sections
        num_moments,max_cell_dof_count);

  return sweep_chunk;
}