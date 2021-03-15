#include "volmesher_extruder.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/** Creates actual z-levels for the input layer specification.*/
void chi_mesh::VolumeMesherExtruder::SetupLayers(int default_layer_count)
{
  //================================================== Create default layers if no
  //                                                   input layers are provided
  if (input_layers.empty())
  {
    chi_log.Log(LOG_0WARNING)
      << "VolumeMesherExtruder: No extrusion layers have been specified. "
      << "A default single layer will be used with height 1.0 and a single "
      << "subdivision.";
    double dz = 1.0/default_layer_count;
    for (int k=0; k<=default_layer_count; k++)
    {
      vertex_layers.push_back(k*dz);
    }
  }
  else
  {
    double last_z=0.0;
    vertex_layers.push_back(last_z);

    for (const auto& input_layer : input_layers)
    {
      double dz = input_layer.height/input_layer.sub_divisions;

      for (int k=0; k<input_layer.sub_divisions; k++)
      {
        last_z += dz;
        vertex_layers.push_back(last_z);
      }
    }
  }

  chi_log.Log(LOG_0)
    << "VolumeMesherExtruder: Total number of cell layers is "
    << vertex_layers.size()-1;
}

