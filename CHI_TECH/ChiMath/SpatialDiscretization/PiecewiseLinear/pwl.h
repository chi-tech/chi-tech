#ifndef _chi_discretization_pwl_h
#define _chi_discretization_pwl_h

#include"ChiMath/SpatialDiscretization/spatial_discretization.h"
#include"../../../ChiMesh/Region/chi_region.h"
#include "CellViews/pwl_cellbase.h"
#include "../../Quadratures/quadrature_triangle.h"
#include "../../Quadratures/quadrature_tetrahedron.h"





//######################################################### Class def
/**Generalization of the Galerkin Finite Element Method
 * with piecewise linear basis functions
 * for use by either a Continues Finite Element Method (CFEM)
 * or a Discontinuous Finite Element Method (DFEM). */
class SpatialDiscretization_PWL : public SpatialDiscretization
{
public:
  std::vector<CellFEView*> cell_fe_views;

private:
  std::vector<bool>        cell_view_added_flags;
  bool mapping_initialized;
public:
  chi_math::QuadratureTriangle*    tri_quad_deg5;
  chi_math::QuadratureTriangle*    tri_quad_deg3_surf;
  chi_math::QuadratureTetrahedron* tet_quad_deg1;
  chi_math::QuadratureTetrahedron* tet_quad_deg3;
  chi_math::QuadratureTetrahedron* tet_quad_deg3_surface;

private:
  std::vector<chi_mesh::Cell*> neighbor_cells;
  std::vector<CellFEView*> neighbor_cell_fe_views;

public:
  //00
  SpatialDiscretization_PWL(int dim=0);
  //01
  void AddViewOfLocalContinuum(chi_mesh::MeshContinuum* grid) override;
  void AddViewOfNeighborContinuums(chi_mesh::MeshContinuum* grid);
  //02
  std::pair<int,int> OrderNodesCFEM(chi_mesh::MeshContinuum* grid);
  CellFEView* MapFeViewL(int cell_local_index);
  int         MapCFEMDOF(int vertex_id);

  //03
  void BuildCFEMSparsityPattern(chi_mesh::MeshContinuum* grid,
                                std::vector<int>& nodal_nnz_in_diag,
                                std::vector<int>& nodal_nnz_off_diag,
                                const std::pair<int,int>& domain_ownership);

  //04
  std::pair<int,int> OrderNodesDFEM(chi_mesh::MeshContinuum* grid);
  int MapDFEMDOF(chi_mesh::Cell* cell, int dof,
                 int component=0,
                 int component_block_offset=1);
  int MapDFEMDOFLocal(chi_mesh::Cell* cell, int dof,
                      int component=0,
                      int component_block_offset=1);

  //05
  void BuildDFEMSparsityPattern(chi_mesh::MeshContinuum* grid,
                                std::vector<int>& nodal_nnz_in_diag,
                                std::vector<int>& nodal_nnz_off_diag,
                                const std::pair<int,int>& domain_ownership);
  chi_mesh::Cell* MapNeighborCell(int cell_glob_index);
  CellFEView* MapNeighborCellFeView(int cell_glob_index);
};

#endif