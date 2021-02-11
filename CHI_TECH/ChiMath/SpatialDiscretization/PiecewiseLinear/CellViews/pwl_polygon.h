#ifndef PWL_POLYGON_VALUES_H
#define PWL_POLYGON_VALUES_H

#include "../pwl.h"
#include <vector>
#include <ChiMesh/Cell/cell_polygon.h>

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
  chi_mesh::Matrix3x3 J;
  chi_mesh::Matrix3x3 Jinv;
  chi_mesh::Matrix3x3 JTinv;
  std::vector<FEqp_data2d> qp_data;
};

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
  std::vector<FEside_data2d> sides;
  chi_math::QuadratureTriangle&      default_volume_quadrature;
  chi_math::QuadratureGaussLegendre& default_surface_quadrature;
public:
  int      num_of_subtris;
  double   beta;
  chi_mesh::Vertex vc;
  std::vector<double> detJ;
  std::vector<int*> node_to_side_map;



public:
  std::vector<chi_mesh::Vector3>                 IntV_gradshapeI;


private:
  chi_mesh::MeshContinuumPtr grid;

  bool precomputed;
  
public:
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constructor
  PolygonPWLFEValues(chi_mesh::CellPolygon* poly_cell,
                     chi_mesh::MeshContinuumPtr vol_continuum,
                     SpatialDiscretization_PWL *discretization);

  double ShapeValue(int i, const chi_mesh::Vector3& xyz) override;
  chi_mesh::Vector3 GradShapeValue(int i, const chi_mesh::Vector3& xyz) override;


  void ShapeValues(const chi_mesh::Vector3& xyz,
                   std::vector<double>& shape_values) override;

  void GradShapeValues(const chi_mesh::Vector3& xyz,
                       std::vector<chi_mesh::Vector3>& gradshape_values) override;


  //############################################### Precomputation cell matrices
  double SideShape(int side, int i, int qpoint_index, bool on_surface = false);
  double SideGradShape_x(int side, int i);
  double SideGradShape_y(int side, int i);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Determinant J
  double DetJ(int s, int qpoint_index, bool on_surface=false)
  {
    if (!on_surface)
      return sides[s].detJ;
    else
      return sides[s].detJ_surf;
  }


public:
  void PreComputeValues() override;

};

#endif