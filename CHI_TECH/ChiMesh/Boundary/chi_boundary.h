#ifndef _chi_boundary_h
#define _chi_boundary_h
#include"../SurfaceMesh/chi_surfacemesh.h"
#include "../MeshContinuum/chi_meshcontinuum.h"

//######################################################### Class definition

class chi_mesh::Boundary
{
public:
  chi_mesh::MeshContinuum* initial_mesh_continuum;
  std::vector<chi_mesh::MeshContinuum*> mesh_continua;


public:
  //00
  Boundary() : initial_mesh_continuum(new chi_mesh::MeshContinuum()) {};

};


#endif