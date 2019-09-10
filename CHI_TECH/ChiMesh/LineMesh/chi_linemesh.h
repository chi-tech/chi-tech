#ifndef _chi_linemesh_h
#define _chi_linemesh_h

#include<stdio.h>
#include <vector>

#include"../chi_mesh.h"

//###################################################################
/**Generalized mesh fulfilling many functions. Primarily this
 * serves as the platform for 1D meshes.*/
class chi_mesh::LineMesh
{
public:
  std::vector<chi_mesh::Vertex> vertices;

};


#endif