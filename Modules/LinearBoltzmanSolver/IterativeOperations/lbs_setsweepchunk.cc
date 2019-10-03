#include "../lbs_linear_boltzman_solver.h"
#include "../SweepChunks/lbs_sweepchunk_pwl_slab.h"
#include "../SweepChunks/lbs_sweepchunk_pwl_polygon.h"
#include "../SweepChunks/lbs_sweepchunk_pwl_polyhedron.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/chi_volumemesher.h>
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include <ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h>
#include <ChiMesh/VolumeMesher/Predefined2D/volmesher_predefined2d.h>

typedef chi_mesh::SweepManagement::SweepChunk SweepChunk;

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
  if      (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
  {
    sweep_chunk =
      new LBSSweepChunkPWLSlab(
        grid,                                    //Spatial grid of cells
        (SpatialDiscretization_PWL*)discretization, //Spatial discretization
        &cell_transport_views,                   //Cell transport views
        &phi_new_local,                          //Destination phi
        &q_moments_local,                        //Source moments
        groupset,                                //Reference groupset
        &material_xs,                            //Material cross-sections
        num_moments,max_cell_dof_count);
  }
  else if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherPredefined2D))
  {
    sweep_chunk =
      new LBSSweepChunkPWLPolygon(
        grid,                                    //Spatial grid of cells
        (SpatialDiscretization_PWL*)discretization, //Spatial discretization
        &cell_transport_views,                   //Cell transport views
        &phi_new_local,                          //Destination phi
        &q_moments_local,                        //Source moments
        groupset,                                //Reference groupset
        &material_xs,                            //Material cross-sections
        num_moments,max_cell_dof_count);
  }
  else if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherExtruder))
  {
    sweep_chunk =
      new LBSSweepChunkPWLPolyhedron(
        grid,                                    //Spatial grid of cells
        (SpatialDiscretization_PWL*)discretization, //Spatial discretization
        &cell_transport_views,                   //Cell transport views
        &phi_new_local,                          //Destination phi
        &q_moments_local,                        //Source moments
        groupset,                                //Reference groupset
        &material_xs,                            //Material cross-sections
        num_moments,max_cell_dof_count);
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "LinearBoltzman::Solver::SetSweepChunk, failed. Could not establish "
      << "which sweep chunk to use.";
    exit(EXIT_FAILURE);
  }


  return sweep_chunk;
}