#ifndef _pwl_polyhedron_h
#define _pwl_polyhedron_h

#include "../pwl.h"
#include <vector>
#include "ChiMath/Quadratures/quadrature.h"
#include "ChiMesh/Cell/cell_polyhedron.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#define ON_SURFACE true

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
  chi_mesh::Vector3         side_centroid;
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
};


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
class PolyhedronFEView : public CellFEView
{
private:
  std::vector<double>            face_betaf;     ///< Face Beta-factor.
  double                         alphac;         ///< Cell alpha-factor.
public:
  std::vector<FEface_data>       face_data;      ///< Holds determinants and data tet-by-tet.
  std::vector<FEnodeMap>         node_side_maps; ///< Maps nodes to side tets.

private:
  std::vector<chi_math::QuadratureTetrahedron*> quadratures; ///< Quadratures used by this method.
  chi_mesh::MeshContinuum*       grid;                       ///< Pointer to the reference grid.

  bool                   precomputed = false;   ///< Are the integrals computed.


public:
  PolyhedronFEView(chi_mesh::CellPolyhedron* polyh_cell,
                   chi_mesh::MeshContinuum* vol_continuum,
                   SpatialDiscretization_PWL* discretization= nullptr);


  //################################################## Define standard
  //                                                   tetrahedron shape
  //                                                   functions
private:
  double TetShape(int index, int qpoint_index, bool on_surface = false);
  static double TetGradShape_x(int index);
  static double TetGradShape_y(int index);
  static double TetGradShape_z(int index);

  //################################################## Shape functions per side
private:
  /**Determinant evaluated at quadrature point*/
  double DetJ(int face_index, int side_index,
              int qpoint_index, bool on_surface=false)
  {
    if (on_surface)
      return (face_data[face_index].sides[side_index].detJ_surf);
    else
      return (face_data[face_index].sides[side_index].detJ);
  }
  /**Shape function evaluation on a tet at a quadrature point*/
  double GetShape(int face, int side, int i, int qp, bool surface = false)
  {
    if (surface)
      return face_data[face].sides[side].qp_data[i].shape_qp_surf[qp];
    else
      return face_data[face].sides[side].qp_data[i].shape_qp[qp];
  }
  /**GradeShape-x function evaluation on a tet at a quadrature point*/
  double GetGradShape_x(int face, int side, int i, int qp)
  {
    return face_data[face].sides[side].qp_data[i].gradshapex_qp[qp];
  }
  /**GradeShape-y function evaluation on a tet at a quadrature point*/
  double GetGradShape_y(int face, int side, int i, int qp)
  {
    return face_data[face].sides[side].qp_data[i].gradshapey_qp[qp];
  }
  /**GradeShape-z function evaluation on a tet at a quadrature point*/
  double GetGradShape_z(int face, int side, int i, int qp)
  {
    return face_data[face].sides[side].qp_data[i].gradshapez_qp[qp];
  }

  //############################################### Actual shape functions
public:
  double ShapeValue(int i, const chi_mesh::Vector3& xyz) override;
  chi_mesh::Vector3 GradShapeValue(int i, const chi_mesh::Vector3& xyz) override;


  void ShapeValues(const chi_mesh::Vector3& xyz,
                   std::vector<double>& shape_values) override;

  void GradShapeValues(const chi_mesh::Vector3& xyz,
                       std::vector<chi_mesh::Vector3>& gradshape_values) override;

  //############################################### Precomputation cell matrices
private:
  double PreShape(int face_index, int side_index,
                  int i, int qpoint_index, bool on_surface = false);

  double PreGradShape_x(int face_index, int side_index,
                        int i);

  double PreGradShape_y(int face_index, int side_index,
                        int i);

  double PreGradShape_z(int face_index, int side_index,
                        int i);

public:
  //####################################################### Precomputing
  void PreCompute();

public:
  void CleanUp()
  {
    for (auto& face : face_data)
      for (auto& side : face.sides)
        side.qp_data = std::move(std::vector<FEqp_data3d>(0));
  }

};

#endif