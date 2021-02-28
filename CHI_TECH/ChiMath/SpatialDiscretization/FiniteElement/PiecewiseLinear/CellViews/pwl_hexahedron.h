#ifndef PWL_HEXAHEDRON_H
#define PWL_HEXAHEDRON_H

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMesh/Cell/cell_polyhedron.h"

#include <vector>

class HexahedronPWLFEValues : public CellMappingFEPWL
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
  const double _1_8th = 1/8.0;
  std::array<double,8> x_i = {0,0,0,0,0,0,0,0},
                       y_i = {0,0,0,0,0,0,0,0},
                       z_i = {0,0,0,0,0,0,0,0};
  static constexpr std::array<double,8>
    xi_i   = {-1, 1, 1, -1, -1, 1, 1, -1},
    eta_i  = {-1, -1, 1, 1, -1, -1, 1, 1},
    zeta_i = {-1, -1, -1, -1, 1, 1, 1, 1};

  static std::vector<MatDbl> qp_df_dg_default;
  static std::vector<MatDbl> qp_df_dg_arbitrary;

  chi_math::QuadratureHexahedron&    default_volume_quadrature;
  chi_math::QuadratureQuadrilateral& default_surface_quadrature;

  chi_math::QuadratureHexahedron&    arbitrary_volume_quadrature;
  chi_math::QuadratureQuadrilateral& arbitrary_surface_quadrature;

public:
  HexahedronPWLFEValues(chi_mesh::CellPolyhedron* polyh_cell,
                        chi_mesh::MeshContinuumPtr ref_grid,
                        chi_math::QuadratureHexahedron&    minumum_volume_quadrature,
                        chi_math::QuadratureQuadrilateral& minumum_surface_quadrature,
                        chi_math::QuadratureHexahedron&    arb_volume_quadrature,
                        chi_math::QuadratureQuadrilateral& arb_surface_quadrature);

  void ComputeUnitIntegrals(
    chi_math::finite_element::UnitIntegralData& ui_data) override;
  void InitializeAllQuadraturePointData(
    chi_math::finite_element::InternalQuadraturePointData& internal_data,
    std::vector<chi_math::finite_element::FaceQuadraturePointData>& faces_qp_data) override;


  //################################################## Define standard
  //                                                   hexahedron linear shape
  //                                                   functions
  //01a_reftet.cc
private:
  double HexShape(unsigned int index, double xi, double eta, double zeta);
  double HexGradShape_x(unsigned int index, double xi, double eta, double zeta);
  double HexGradShape_y(unsigned int index, double xi, double eta, double zeta);
  double HexGradShape_z(unsigned int index, double xi, double eta, double zeta);

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


#endif //PWL_HEXAHEDRON_H