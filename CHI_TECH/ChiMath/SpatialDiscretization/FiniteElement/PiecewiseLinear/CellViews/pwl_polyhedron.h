#ifndef PWL_POLYHEDRON_VALUES_H
#define PWL_POLYHEDRON_VALUES_H

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include <vector>
#include "ChiMath/Quadratures/quadrature.h"
#include "ChiMesh/Cell/cell_polyhedron.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"



//###################################################################
/**Object for handling piecewise linear
 * shape functions on polyhedron shaped 3D cells.
 *
 * This object has a whitepaper associated with it
 * (<a target="_blank"
 * href="../../whitepages/FEM/PWLPolyhedron/PWLPolyhedron.pdf">
 * here</a>)
 *
 * Some notes on indexing:\n
 *  - IntV_shapeI_gradshapeJ, given i and j results in a vector.
 *  - IntS_shapeI_shapeJ, requires f, then cell i, then cell j
 *  - face_dof_mappings, is as follows face_dof_mappings[f][fi] and
 *    returns cell i.
 * */
class PolyhedronPWLFEValues : public CellPWLFEValues
{
private:
  /**For a given side(tet), this structure holds the values of
 * shape functions at each quadrature point.*/
  struct FEqp_data3d
  {
    std::vector<double> shape_qp;
    std::vector<double> shape_qp_surf;
    std::vector<double> gradshapex_qp;
    std::vector<double> gradshapey_qp;
    std::vector<double> gradshapez_qp;
  };
  //Goes into
  /**Stores the data for each side's tetrahedron. */
  struct FEside_data3d
  {
    double                    detJ = 0.0;
    double                    detJ_surf = 0.0;
    std::vector<int>          v_index;
    chi_mesh::Vector3         v0;
    chi_mesh::Matrix3x3       J;
    chi_mesh::Matrix3x3       Jinv;
    chi_mesh::Matrix3x3       JTinv;
    std::vector<FEqp_data3d>  qp_data;
  };
  //Goes into
  /**Stores data for each face.*/
  struct FEface_data
  {
    std::vector<FEside_data3d> sides;
    chi_mesh::Vector3 normal;
  };
  //Goes int face_data


  /**Lowest level of mapping dof i.*/
  struct FEnodeSideMap
  {
    int index = -1;
    bool part_of_face = false;
  };
  //Goes into
  /**Intermediate level of mapping.*/
  struct FEnodeFaceMap
  {
    std::vector<FEnodeSideMap> side_map;
  };
  //Goes into
  /**Node map per face.*/
  struct FEnodeMap
  {
    std::vector<FEnodeFaceMap> face_map;
  };
  //Goes into node_maps
  // node n
  // face f
  // side s
  // node_maps[n]->face_map[f]->side_map[s]

private:
  std::vector<double>            face_betaf;     ///< Face Beta-factor.
  double                         alphac;         ///< Cell alpha-factor.

  std::vector<FEface_data>       face_data;      ///< Holds determinants and data tet-by-tet.
private:
  std::vector<FEnodeMap>         node_side_maps; ///< Maps nodes to side tets.

private:
  chi_math::QuadratureTetrahedron& default_volume_quadrature;
  chi_math::QuadratureTriangle&    default_surface_quadrature;

  chi_math::QuadratureTetrahedron& arbitrary_volume_quadrature;
  chi_math::QuadratureTriangle&    arbitrary_surface_quadrature;

  chi_math::QuadratureTetrahedron* active_volume_quadrature =nullptr;
  chi_math::QuadratureTriangle*    active_surface_quadrature=nullptr;

public:
  //00_constrdestr.cc
  PolyhedronPWLFEValues(chi_mesh::CellPolyhedron* polyh_cell,
                        chi_mesh::MeshContinuumPtr ref_grid,
                        chi_math::QuadratureTetrahedron& minumum_volume_quadrature,
                        chi_math::QuadratureTriangle&    minumum_surface_quadrature,
                        chi_math::QuadratureTetrahedron& arb_volume_quadrature,
                        chi_math::QuadratureTriangle&    arb_surface_quadrature);

  void ComputeUnitIntegrals(
    chi_math::finite_element::UnitIntegralData& ui_data) override;
  void InitializeQuadraturePointData(
    chi_math::finite_element::InternalQuadraturePointData& internal_data,
    std::vector<chi_math::finite_element::FaceQuadraturePointData>& faces_qp_data) override;

  //################################################## Define standard
  //                                                   tetrahedron linear shape
  //                                                   functions
  //01a_reftet.cc
private:
  double TetShape(int index, int qpoint_index, bool on_surface = false);
  static double TetGradShape_x(int index);
  static double TetGradShape_y(int index);
  static double TetGradShape_z(int index);

  //################################################## Shape functions per face-side
  //01b_sidevalues.cc
private:
  double FaceSideShape(int face_index, int side_index, int i,
                       int qpoint_index, bool on_surface = false);
  double FaceSideGradShape_x(int face_index, int side_index, int i);
  double FaceSideGradShape_y(int face_index, int side_index, int i);
  double FaceSideGradShape_z(int face_index, int side_index, int i);

  //############################################### Actual shape functions
  //                                                as function of cartesian
  //                                                coordinates
public:
  double ShapeValue(int i, const chi_mesh::Vector3& xyz) override;
  chi_mesh::Vector3 GradShapeValue(int i, const chi_mesh::Vector3& xyz) override;


  void ShapeValues(const chi_mesh::Vector3& xyz,
                   std::vector<double>& shape_values) override;

  void GradShapeValues(const chi_mesh::Vector3& xyz,
                       std::vector<chi_mesh::Vector3>& gradshape_values) override;

};

#endif
