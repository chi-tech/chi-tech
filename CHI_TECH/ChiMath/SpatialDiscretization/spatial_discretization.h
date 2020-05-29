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

  chi_mesh::MeshContinuum* ref_grid;

public:
  std::vector<int> node_mapping;
  std::vector<int> reverse_node_mapping;
  int              fv_local_block_address = 0;
  int              cfem_local_block_address = 0;
  int              dfem_local_block_address = 0;
  std::vector<int> cell_dfem_block_address;
  std::vector<std::pair<int,int>> neighbor_cell_block_address;

public:
  int              block_size_per_unknown=0;

public:
  //00
  SpatialDiscretization(int dim);

  //01
  virtual void AddViewOfLocalContinuum(chi_mesh::MeshContinuum* grid);

};

#endif