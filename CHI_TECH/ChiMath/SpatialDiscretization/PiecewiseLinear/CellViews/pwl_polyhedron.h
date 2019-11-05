#ifndef _pwl_polyhedron_h
#define _pwl_polyhedron_h

#include "../pwl.h"
#include <vector>
#include "../../../Quadratures/quadrature.h"
#include "../../../../ChiMesh/Cell/cell_polyhedron.h"
#include "ChiMesh/Cell/cell_polyhedronv2.h"
#include "../../../../ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#define ON_SURFACE true

/**For a given side(tet), this structure holds the values of
 * shape functions at each quadrature point.*/
struct FEqp_data3d
{
  std::vector<float> shape_qp;
  std::vector<float> shape_qp_surf;
  std::vector<float> gradshapex_qp;
  std::vector<float> gradshapey_qp;
  std::vector<float> gradshapez_qp;
};
//Goes into
/**Stores the data for each side's tetrahedron. */
struct FEside_data3d
{
  float detJ;
  float detJ_surf;
  int*    v_index;
  chi_mesh::Vector sc;
  chi_mesh::Matrix3x3 J;
  chi_mesh::Matrix3x3 Jinv;
  chi_mesh::Matrix3x3 JTinv;
  std::vector<FEqp_data3d*> qp_data;
};
//Goes into
/**Stores data for each face.*/
struct FEface_data
{
  std::vector<FEside_data3d*> sides;
  chi_mesh::Vector vfc;
};


/**Lowest level of mapping dof i.*/
struct FEnodeSideMap
{
  int index;
  bool part_of_face;
};
//Goes into
/**Intermediate level of mapping.*/
struct FEnodeFaceMap
{
  std::vector<FEnodeSideMap*> side_map;
};
//Goes into
/**Node map per face.*/
struct FEnodeMap
{
  std::vector<FEnodeFaceMap*> face_map;
};
//Goes into node_maps
// node n
// face f
// side s
// node_maps[n]->face_map[f]->side_map[s]


/**FaceDOFmapping*/
//struct FEFaceDOFMapping
//{
//  std::vector<int> cell_dof;
//};

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
  chi_mesh::Vertex                vcc;                   ///< Cell centroid
  std::vector<FEface_data*>       faces;
  std::vector<double>              face_betaf;
public:
  std::vector<FEnodeMap*>         node_maps;

private:
  std::vector<std::vector<std::vector<double>>> IntSi_shapeI_shapeJ;
  std::vector<std::vector<std::vector<chi_mesh::Vector>>> IntSi_shapeI_gradshapeJ;

private:
  std::vector<chi_math::QuadratureTetrahedron*> quadratures;
  chi_mesh::MeshContinuum* grid;

  double alphac;
  bool precomputed;


public:
  PolyhedronFEView(chi_mesh::CellPolyhedron* polyh_cell,
                   chi_mesh::MeshContinuum* vol_continuum,
                   SpatialDiscretization_PWL* discretization= nullptr);
  PolyhedronFEView(chi_mesh::CellPolyhedronV2* polyh_cell,
                   chi_mesh::MeshContinuum* vol_continuum,
                   SpatialDiscretization_PWL* discretization= nullptr);


  //################################################## Define standard
  //                                                   tetrahedron shape
  //                                                   functions
private:
  double TetShape(int index, int qpoint_index, bool on_surface = false);
  double TetGradShape_x(int index, int qpoint_index);
  double TetGradShape_y(int index, int qpoint_index);
  double TetGradShape_z(int index, int qpoint_index);

  //################################################## Shape functions per side
private:
  /**Determinant evaluated at quadrature point*/
  double DetJ(int face_index, int side_index,
              int qpoint_index, bool on_surface=false)
  {
    if (on_surface)
      return (faces[face_index]->sides[side_index]->detJ_surf);
    else
      return (faces[face_index]->sides[side_index]->detJ);
  }
  double GetShape(int face, int side, int i, int qp, bool surface = false)
  {
    if (surface)
      return faces[face]->sides[side]->qp_data[i]->shape_qp_surf[qp];
    else
      return faces[face]->sides[side]->qp_data[i]->shape_qp[qp];
  }

  double GetGradShape_x(int face, int side, int i, int qp)
  {
    return faces[face]->sides[side]->qp_data[i]->gradshapex_qp[qp];
  }

  double GetGradShape_y(int face, int side, int i, int qp)
  {
    return faces[face]->sides[side]->qp_data[i]->gradshapey_qp[qp];
  }

  double GetGradShape_z(int face, int side, int i, int qp)
  {
    return faces[face]->sides[side]->qp_data[i]->gradshapez_qp[qp];
  }

  //############################################### Actual shape functions
public:
  double           Shape_xyz(int i, chi_mesh::Vector xyz);
  chi_mesh::Vector GradShape_xyz(int i, chi_mesh::Vector xyz);

  //############################################### Precomputation cell matrices
private:
  double PreShape(int face_index, int side_index,
                  int i, int qpoint_index, bool on_surface = false);

  double PreGradShape_x(int face_index, int side_index,
                        int i, int qpoint_index);

  double PreGradShape_y(int face_index, int side_index,
                        int i, int qpoint_index);

  double PreGradShape_z(int face_index, int side_index,
                        int i, int qpoint_index);

public:
  //####################################################### Precomputing
  void PreCompute();

public:
  void CleanUp()
  {
//    face_betaf.clear(); face_betaf.shrink_to_fit();
//    for (int f=(faces.size()-1); f>=0; f--)
//    {
//      for (int s=(faces[f]->sides.size()-1); s>=0; s--)
//      {
//        //delete [] faces[f]->sides[s]->m;
//        delete [] faces[f]->sides[s]->v_index;
//
//        for (int qp=(faces[f]->sides[s]->qp_data.size()-1); qp>=0; qp--)
//        {
//          FEqp_data* cur_qp = faces[f]->sides[s]->qp_data[qp];
//          cur_qp->varphi_qp.clear();
//          cur_qp->gradvarphix_qp.clear();
//          cur_qp->gradvarphiy_qp.clear();
//          cur_qp->gradvarphiz_qp.clear();
//
//          cur_qp->varphi_qp.shrink_to_fit();
//          cur_qp->gradvarphix_qp.shrink_to_fit();
//          cur_qp->gradvarphiy_qp.shrink_to_fit();
//          cur_qp->gradvarphiz_qp.shrink_to_fit();
//
//          faces[f]->sides[s]->qp_data.erase(
//            faces[f]->sides[s]->qp_data.begin()+qp);
//          delete cur_qp;
//        }
//        FEside_data* cur_side = faces[f]->sides[s];
//        delete cur_side;
//      }
//      FEface_data* cur_face = faces[f];
//      delete cur_face;
//    }
//    faces.clear();
//    faces.shrink_to_fit();
//
//    for (int n=0; n<node_maps.size(); n++)
//    {
//      FEnodeMap* cur_node = node_maps[n];
//      for (int f=0; f<cur_node->face_map.size(); f++)
//      {
//        FEnodeFaceMap* cur_face = cur_node->face_map[f];
//        for (int s=0; s<node_maps[n]->face_map[f]->side_map.size(); s++)
//        {
//          FEnodeSideMap* cur_side = cur_face->side_map[s];
//          delete [] cur_side;
//        }
//        delete cur_face;
//      }
//      delete cur_node;
//    }
//    node_maps.clear();
//    node_maps.shrink_to_fit();
  }

};

#endif