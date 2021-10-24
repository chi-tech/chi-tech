#ifndef PWL_POLYGON_VALUES_H
#define PWL_POLYGON_VALUES_H

#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"
#include "ChiMath/Quadratures/quadrature_line.h"
#include "ChiMath/Quadratures/quadrature_triangle.h"
#include "ChiMesh/Cell/cell.h"



//###################################################################
/** Object for handling polygon shaped 2D cells.
 *
 * This object has a whitepaper associated with it
 * (<a target="_blank"
 * href="../../whitepages/FEM/PWLPolygon/PWLPolygon.pdf">
 * here</a>)
 *
 * Notes on indexing:\n
 * - IntS_shapeI_shapeJ is indexed as [f][i][j]
 * - IntS_shapeI is indexed as [i][f]
 * - IntS_shapeI_gradshapeJ is indexed as [f][i][j]
 * - node_to_side_map is indexed as [i][f]
 * - edge_dof_mappings, is indexed as [f][fi] and
 *    returns cell dof i.*/
class PolygonMappingFE_PWL : public CellMappingFE_PWL
{
private:
  struct FEside_data2d
  {
    double detJ;
    double detJ_surf;
    std::array<int,2> v_index;
    chi_mesh::Vector3   v0;
    chi_mesh::Matrix3x3 J;
    chi_mesh::Matrix3x3 Jinv;
    chi_mesh::Matrix3x3 JTinv;
    chi_mesh::Vector3 normal;
  };
  //Goes into sides

  std::vector<FEside_data2d> sides;
  const chi_math::QuadratureTriangle& volume_quadrature;
  const chi_math::QuadratureLine&     surface_quadrature;

private:
  int      num_of_subtris;
  double   beta;
  chi_mesh::Vertex vc;
  std::vector<std::vector<int>> node_to_side_map;



public:
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constructor
  PolygonMappingFE_PWL(const chi_mesh::Cell& poly_cell,
                       const chi_mesh::MeshContinuumPtr& ref_grid,
                       const chi_math::QuadratureTriangle& volume_quadrature,
                       const chi_math::QuadratureLine&     surface_quadrature);

  void InitializeVolumeQuadraturePointData(
    chi_math::finite_element::InternalQuadraturePointData& internal_data) const override;

  void InitializeFaceQuadraturePointData(
    unsigned int face,
    chi_math::finite_element::FaceQuadraturePointData& faces_qp_data) const override;

  //################################################## Define standard
  //                                                   triangle linear shape
  //                                                   functions
  static
  double TriShape(int index,
                  const chi_mesh::Vector3& qpoint,
                  bool on_surface = false);

  //############################################### Shape functions per side
  double SideShape(unsigned int side,
                   unsigned int i,
                   const chi_mesh::Vector3& qpoint,
                   bool on_surface = false) const;
  double SideGradShape_x(unsigned int side, int i) const;
  double SideGradShape_y(unsigned int side, int i) const;

  double ShapeValue(int i, const chi_mesh::Vector3& xyz) override;
  chi_mesh::Vector3 GradShapeValue(int i, const chi_mesh::Vector3& xyz) override;


  void ShapeValues(const chi_mesh::Vector3& xyz,
                   std::vector<double>& shape_values) override;

  void GradShapeValues(const chi_mesh::Vector3& xyz,
                       std::vector<chi_mesh::Vector3>& gradshape_values) override;
};

#endif