#ifndef SPATIAL_DISCRETIZATION_PWLC_H
#define SPATIAL_DISCRETIZATION_PWLC_H

#include "CellViews/pwl_cellbase.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/Quadratures/quadrature_gausslegendre.h"
#include "ChiMath/Quadratures/quadrature_triangle.h"
#include "ChiMath/Quadratures/quadrature_tetrahedron.h"

//######################################################### Class def
/**Generalization of the Galerkin Finite Element Method
 * with piecewise linear basis functions
 * for use by either a Continues Finite Element Method (CFEM)
 * or a Discontinuous Finite Element Method (DFEM). */
class SpatialDiscretization_PWLC : public SpatialDiscretization
{
public:
  std::vector<CellPWLFEValues*> cell_fe_views;

private:
  std::vector<bool>        cell_view_added_flags;
  bool                     mapping_initialized;
public:
  chi_math::QuadratureGaussLegendre line_quad_order_second;
  chi_math::QuadratureTriangle      tri_quad_order_second;
  chi_math::QuadratureTetrahedron   tet_quad_order_second;

  std::map<int,int> node_mapping;

  int local_block_address = 0;
//  std::vector<int> cell_local_block_address;
//  std::vector<std::pair<int,int>> neighbor_cell_block_address;

  std::vector<int> locJ_block_address;
  std::vector<int> locJ_block_size;

  unsigned int local_base_block_size=0;
  unsigned int globl_base_block_size=0;

private:
  std::vector<chi_mesh::Cell*> neighbor_cells;
  std::vector<CellPWLFEValues*> neighbor_cell_fe_views;

private:
  //00
  explicit
  SpatialDiscretization_PWLC(chi_mesh::MeshContinuumPtr in_grid);

public:
  //prevent anything else other than a shared pointer
  static
  std::shared_ptr<SpatialDiscretization_PWLC>
  New(chi_mesh::MeshContinuumPtr in_grid)
  { if (in_grid == nullptr) throw std::invalid_argument(
      "Null supplied as grid to SpatialDiscretization_PWLC.");
    return std::shared_ptr<SpatialDiscretization_PWLC>(
      new SpatialDiscretization_PWLC(in_grid));}

  //01
  void PreComputeCellSDValues(chi_mesh::MeshContinuumPtr grid) override;
  void PreComputeNeighborCellSDValues(chi_mesh::MeshContinuumPtr grid);
  CellPWLFEValues* MapFeViewL(int cell_local_index);

private:
  //02
  void OrderNodes(chi_mesh::MeshContinuumPtr grid);

public:
  //03
  void BuildSparsityPattern(chi_mesh::MeshContinuumPtr grid,
                            std::vector<int>& nodal_nnz_in_diag,
                            std::vector<int>& nodal_nnz_off_diag,
                            chi_math::UnknownManager& unknown_manager) override;
//  chi_mesh::Cell* MapNeighborCell(int cell_glob_index);
//  CellPWLFEValues* MapNeighborCellFeView(int cell_glob_index);

  //04 Mappings
  int MapDOF(int vertex_id);



  int MapDOF(int vertex_id,
             chi_math::UnknownManager& unknown_manager,
             unsigned int unknown_id,
             unsigned int component= 0);

  //05
  size_t GetNumLocalDOFs(chi_mesh::MeshContinuumPtr grid,
                         chi_math::UnknownManager& unknown_manager) override;
  size_t GetNumGlobalDOFs(chi_mesh::MeshContinuumPtr grid,
                          chi_math::UnknownManager& unknown_manager) override;
//  unsigned int GetNumGhostDOFs(chi_mesh::MeshContinuumPtr grid,
//                               chi_math::UnknownManager* unknown_manager);
//
//  std::vector<int> GetGhostDOFIndices(chi_mesh::MeshContinuumPtr grid,
//                                      chi_math::UnknownManager* unknown_manager,
//                                      unsigned int unknown_id=0);

  void LocalizePETScVector(Vec petsc_vector,
                           std::vector<double>& local_vector,
                           chi_math::UnknownManager& unknown_manager)
                           override;
};

#endif //SPATIAL_DISCRETIZATION_PWLC_H