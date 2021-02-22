#ifndef PWL_POLYGON_VALUES_H
#define PWL_POLYGON_VALUES_H

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include <vector>
#include <ChiMesh/Cell/cell_polygon.h>



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
class PolygonPWLFEValues : public CellPWLFEValues
{
private:
  /**For a given side(triangle), this structure holds the values of
 * shape functions at each quadrature point.*/
  struct FEqp_data2d
  {
    std::vector<double> shape_qp;
    std::vector<double> shape_qp_surf;
    std::vector<double> gradshapex_qp;
    std::vector<double> gradshapey_qp;
  };
  //Goes into
  struct FEside_data2d
  {
    double detJ;
    double detJ_surf;
    std::array<int,2> v_index;
    chi_mesh::Vector3   v0;
    chi_mesh::Matrix3x3 J;
    chi_mesh::Matrix3x3 Jinv;
    chi_mesh::Matrix3x3 JTinv;
    std::vector<FEqp_data2d> qp_data;
    chi_mesh::Vector3 normal;
  };
  //Goes into sides

  std::vector<FEside_data2d> sides;
  chi_math::QuadratureTriangle&      default_volume_quadrature;
  chi_math::QuadratureGaussLegendre& default_surface_quadrature;

  chi_math::QuadratureTriangle&      arbitrary_volume_quadrature;
  chi_math::QuadratureGaussLegendre& arbitrary_surface_quadrature;

  chi_math::QuadratureTriangle*      active_volume_quadrature = nullptr;
  chi_math::QuadratureGaussLegendre* active_surface_quadrature= nullptr;
private:
  int      num_of_subtris;
  double   beta;
  chi_mesh::Vertex vc;
  std::vector<double> detJ;
  std::vector<std::vector<int>> node_to_side_map;



//public:
//  std::vector<chi_mesh::Vector3>                 IntV_gradshapeI;
  
public:
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constructor
  PolygonPWLFEValues(chi_mesh::CellPolygon* poly_cell,
                     chi_mesh::MeshContinuumPtr ref_grid,
                     chi_math::QuadratureTriangle&      minumum_volume_quadrature,
                     chi_math::QuadratureGaussLegendre& minumum_surface_quadrature,
                     chi_math::QuadratureTriangle&      arb_volume_quadrature,
                     chi_math::QuadratureGaussLegendre& arb_surface_quadrature);

  void ComputeUnitIntegrals(
    chi_math::finite_element::UnitIntegralData& ui_data) override;
  void InitializeQuadraturePointData(
    chi_math::finite_element::InternalQuadraturePointData& internal_data,
    std::vector<chi_math::finite_element::FaceQuadraturePointData>& faces_qp_data) override;

  //################################################## Define standard
  //                                                   triangle linear shape
  //                                                   functions
  double TriShape(int index, int qpoint_index, bool on_surface = false);

  //############################################### Shape functions per side
  double SideShape(int side, int i, int qpoint_index, bool on_surface = false);
  double SideGradShape_x(int side, int i);
  double SideGradShape_y(int side, int i);

  double ShapeValue(int i, const chi_mesh::Vector3& xyz) override;
  chi_mesh::Vector3 GradShapeValue(int i, const chi_mesh::Vector3& xyz) override;


  void ShapeValues(const chi_mesh::Vector3& xyz,
                   std::vector<double>& shape_values) override;

  void GradShapeValues(const chi_mesh::Vector3& xyz,
                       std::vector<chi_mesh::Vector3>& gradshape_values) override;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Determinant J
  double DetJ(int s, int qpoint_index, bool on_surface=false)
  {
    if (!on_surface)
      return sides[s].detJ;
    else
      return sides[s].detJ_surf;
  }

};

#endif