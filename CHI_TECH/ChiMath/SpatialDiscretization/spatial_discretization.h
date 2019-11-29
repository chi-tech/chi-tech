#ifndef _chi_discretization_h
#define _chi_discretization_h

#include "../../ChiMesh/chi_mesh.h"
#include "../Quadratures/quadrature.h"

#define DEG1         0
#define DEG3         1
#define DEG3_SURFACE 2


class SpatialDiscretization
{
public:
  int dim;

public:
  std::vector<int> node_mapping;
  std::vector<int> reverse_node_mapping;

public:
  //00
  SpatialDiscretization(int dim);

  //01
  virtual void AddViewOfLocalContinuum(
    chi_mesh::MeshContinuum* vol_continuum,
    int num_cells,
    int* cell_indices);

  virtual void AddViewOfLocalContinuum(chi_mesh::MeshContinuum* vol_continuum);

};

#endif