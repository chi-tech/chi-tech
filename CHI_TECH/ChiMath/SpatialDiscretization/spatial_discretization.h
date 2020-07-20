#ifndef _chi_discretization_h
#define _chi_discretization_h

#include "ChiMesh/chi_mesh.h"
#include "../Quadratures/quadrature.h"
#include "ChiMath/chi_math.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <petscksp.h>

#define DEG1         0
#define DEG3         1
#define DEG3_SURFACE 2


class SpatialDiscretization
{
public:
  int dim;
  const chi_math::SpatialDiscretizationType type;

  chi_mesh::MeshContinuum* ref_grid;

public:
  std::vector<int> node_mapping;
  std::vector<int> reverse_node_mapping;
  int              fv_local_block_address = 0;
  int              cfem_local_block_address = 0;
  int              dfem_local_block_address = 0;
  std::vector<int> cell_dfem_block_address;
  std::vector<std::pair<int,int>> neighbor_cell_block_address;

protected:
  //Pretending that there is only one unknown
  unsigned int local_base_block_size=0;
  unsigned int globl_base_block_size=0;

  std::vector<int> locJ_block_address;
  std::vector<int> locJ_block_size;

public:
  int              block_size_per_unknown=0;

public:
  //00
  SpatialDiscretization(int dim,
                        chi_math::SpatialDiscretizationType in_type =
                          chi_math::SpatialDiscretizationType::UNDEFINED);

  //01
  virtual void AddViewOfLocalContinuum(chi_mesh::MeshContinuum* grid);

  //02
  /**Develops a localized view of a petsc vector.
   * Each spatial discretization has a specialization of this
   * method.*/
  virtual void LocalizePETScVector(Vec petsc_vector,
                                   std::vector<double>& local_vector,
                                   chi_math::UnknownManager* unknown_manager)
  {}

};

#endif