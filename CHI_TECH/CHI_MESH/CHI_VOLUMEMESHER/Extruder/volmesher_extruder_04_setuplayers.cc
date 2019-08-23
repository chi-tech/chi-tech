#include "volmesher_extruder.h"

#include <chi_log.h>
extern CHI_LOG chi_log;

//###################################################################
/** Creates actual z-levels for the input layer specification.*/
void chi_mesh::VolumeMesherExtruder::SetupLayers(int default_layer_count)
{
  //================================================== Create default layers if no
  //                                                   input layers are provided
  if (input_layers.size()==0)
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

    for (int ell=0; ell<input_layers.size(); ell++)
    {
      double dz = input_layers[ell]->height/input_layers[ell]->sub_divisions;

      for (int k=0; k<input_layers[ell]->sub_divisions; k++)
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