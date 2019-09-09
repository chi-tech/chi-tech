#ifndef _chi_boundary_h
#define _chi_boundary_h
#include"../CHI_SURFACEMESH/chi_surfacemesh.h"
#include "../CHI_MESHCONTINUUM/chi_meshcontinuum.h"

//######################################################### Class definition

class chi_mesh::Boundary
{
public:
  chi_mesh::MeshContinuum initial_mesh_continuum;
  std::vector<chi_mesh::MeshContinuum*> mesh_continua;


public:
  //00
  Boundary();

};


#endif