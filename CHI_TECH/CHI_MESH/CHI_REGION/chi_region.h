#ifndef _chi_region_h
#define _chi_region_h

#include "../CHI_MESHCONTINUUM/chi_meshcontinuum.h"



//######################################################### Class definition
class chi_mesh::Region
{
public:
  std::vector<chi_mesh::Boundary*> boundaries;
  std::vector<chi_mesh::MeshContinuum*> volume_mesh_continua;
};

class chi_mesh::EmptyRegion : public chi_mesh::Region
{

};

#endif