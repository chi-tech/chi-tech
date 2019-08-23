#ifndef _chi_discretization_h
#define _chi_discretization_h

#include "../../CHI_MESH/chi_mesh.h"
#include "../Quadratures/quadrature.h"

#define DEG1         0
#define DEG3         1
#define DEG3_SURFACE 2


class CHI_DISCRETIZATION
{
public:
  int dim;

public:
  //00
  CHI_DISCRETIZATION(int dim);

  //01
  virtual void AddViewOfLocalContinuum(
    chi_mesh::MeshContinuum* vol_continuum,
    int num_cells,
    int* cell_indices);

};

#endif