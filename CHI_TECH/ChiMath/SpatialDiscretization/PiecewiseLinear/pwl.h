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
  std::vector<int>         cell_fe_views_mapping;
  bool mapping_initialized;
  chi_math::QuadratureTriangle*    tri_quad_deg5;
  chi_math::QuadratureTriangle*    tri_quad_deg3_surf;
  chi_math::QuadratureTetrahedron* tet_quad_deg1;
  chi_math::QuadratureTetrahedron* tet_quad_deg3;
  chi_math::QuadratureTetrahedron* tet_quad_deg3_surface;

public:
  //00
  SpatialDiscretization_PWL(int dim=0);
  //01
  void AddViewOfLocalContinuum(
    chi_mesh::MeshContinuum* vol_continuum,
    int num_cells,
    int* cell_indices);
  void AddViewOfLocalContinuum(chi_mesh::MeshContinuum* vol_continuum) override;
  //02
  std::pair<int,int> OrderNodesCFEM(chi_mesh::MeshContinuum* grid);
  CellFEView* MapFeView(int cell_glob_index);
  int         MapNode(int vertex_id);

  //03
  void BuildCFEMSparsityPattern(chi_mesh::MeshContinuum* grid,
                                std::vector<int>& nodal_bndry_ids,
                                std::vector<int>& nodal_nnz_in_diag,
                                std::vector<int>& nodal_nnz_off_diag,
                                const std::pair<int,int>& domain_ownership);
};

#endif