#ifndef _chi_discretization_pwl_h
#define _chi_discretization_pwl_h

#include"CHI_MATH/SpatialDiscretization/spatial_discretization.h"
#include"../../../CHI_MESH/CHI_REGION/chi_region.h"
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
  CHI_QUADRATURE_TRIANGLE*    tri_quad_deg5;
  CHI_QUADRATURE_TRIANGLE*    tri_quad_deg3_surf;
  CHI_QUADRATURE_TETRAHEDRON* tet_quad_deg1;
  CHI_QUADRATURE_TETRAHEDRON* tet_quad_deg3;
  CHI_QUADRATURE_TETRAHEDRON* tet_quad_deg3_surface;

public:
  //00
  SpatialDiscretization_PWL(int dim=0);
  //01
  void AddViewOfLocalContinuum(
    chi_mesh::MeshContinuum* vol_continuum,
    int num_cells,
    int* cell_indices);
  CellFEView* MapFeView(int cell_glob_index);
};

#endif