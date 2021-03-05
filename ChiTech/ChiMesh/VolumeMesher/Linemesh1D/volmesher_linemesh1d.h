#ifndef _volmesher_linemesh1D_h
#define _volmesher_linemesh1D_h

#include "../chi_volumemesher.h"

//###################################################################
/**A simple mesher to convert line meshes to slab meshes.*/
class chi_mesh::VolumeMesherLinemesh1D : public chi_mesh::VolumeMesher
{
public:
  int num_slab_cells;
public:
  VolumeMesherLinemesh1D()
  {
    num_slab_cells = 0;
  }
  void Execute() override;
};

#endif