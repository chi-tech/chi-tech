#ifndef _volmesher_predefined2d_h
#define _volmesher_predefined2d_h

#include "../chi_volumemesher.h"
//#include "../../MeshContinuum/chi_meshcontinuum.h"

class chi_mesh::VolumeMesherPredefined2D : public chi_mesh::VolumeMesher
{
public:
  bool always_create_polygons;
public:
  //00
       VolumeMesherPredefined2D();
  //01 Utils

  //02
  void Execute();
  //03

};

#endif