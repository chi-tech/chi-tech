#ifndef _chi_mesh_h
#define _chi_mesh_h

#include<vector>
#include<iostream>



/** # Namespace for all meshing features
 *
 * Meshes in ChiTech follow the concept of Regions. In any given region the
 * boundaries are a collection of either line-meshes (2D) or
 * surface-meshes (3D).
 * */
namespace chi_mesh
{
  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FORWARD DECLARATIONS
  struct Vector3;
  typedef Vector3 Normal;
  typedef Vector3 Vertex;
  typedef Vector3 Node;

  struct Matrix3x3;
  struct TensorRank2Dim3;

  struct Face;     //Triangle
  struct Edge;
  struct EdgeLoop;
  struct PolyFace;

  typedef Edge Line;
  typedef std::vector<Edge>      EdgeList;
  typedef std::vector<EdgeLoop*> EdgeLoopCollection;

  struct CellIndexMap;
  typedef CellIndexMap NodeIndexMap;

  struct CELL_SET;

  struct SweepPlane;
  struct SPDS;

  //=================================== Boundary
  class Boundary;

  //=================================== Cells
  class Cell;
//  class CellSlab;
//  class CellTriangle;
//  class CellPolygon;
//  class CellPolyhedron;

  //=================================== Field function interpolation
  class FieldFunctionInterpolation;
  class FieldFunctionInterpolationSlice;
  class FieldFunctionInterpolationLine;
  class FieldFunctionInterpolationVolume;

  //=================================== Meshes
  class LineMesh;
  class SurfaceMesh;
  class UnpartitionedMesh;
  class MeshContinuum;

  //=================================== Logical Volumes
  class LogicalVolume;
  class SphereLogicalVolume;
  class RPPLogicalVolume;
  class RCCLogicalVolume;
  class SurfaceMeshLogicalVolume;
  class BooleanLogicalVolume;

  //=================================== Mesh handler
  class MeshHandler;

  //=================================== Region
  class Region;
  class EmptyRegion;

  //=================================== Surface Meshers
  class SurfaceMesher;
  class SurfaceMesherPassthrough;
  class SurfaceMesherPredefined;
  class SurfaceMesherDelaunay;
  class SurfaceMesherTriangle;

  struct Interface;

  //==================================== Volume meshers
  class VolumeMesher;
  class VolumeMesherLinemesh1D;
  class VolumeMesherPredefined2D;
  class VolumeMesherExtruder;
  class VolumeMesherPredefined3D;
  class VolumeMesherPredefinedUnpartitioned;




  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ROUTINES

  Boundary*              AssignSurfaceToBoundary(chi_mesh::SurfaceMesh* surface);
  MeshHandler*           GetCurrentHandler();
  size_t                 PushNewHandler();
  MeshHandler*           GetNewHandler();
  EdgeLoopCollection*    SplitEdgeLoopByAngle(EdgeLoop* input,double angle=1);

  //=================================== Domain decompositions
  double ComputeLBF(std::vector<Vector3>& points,
                    std::vector<double>& x_cuts,
                    std::vector<double>& y_cuts);
  void   DecomposeSurfaceMeshPxPy(SurfaceMesh* smesh, int Px, int Py);
  void   SurfaceMeshImprintLine(SurfaceMesh* smesh,
                                Vertex* line_point_0,
                                Vertex* line_point_1,
                                double tolerance);

  //=================================== OrthoMeshes
  void Create1DSlabMesh(std::vector<double>& vertices_1d);

  void Create2DOrthoMesh(std::vector<double>& vertices_1d_x,
                         std::vector<double>& vertices_1d_y);

  void Create3DOrthoMesh(std::vector<double>& vertices_1d_x,
                         std::vector<double>& vertices_1d_y,
                         std::vector<double>& vertices_1d_z);

  void CreateUnpartitioned1DOrthoMesh(std::vector<double>& vertices_1d);

  void CreateUnpartitioned2DOrthoMesh(std::vector<double>& vertices_1d_x,
                                      std::vector<double>& vertices_1d_y);

  void CreateUnpartitioned3DOrthoMesh(std::vector<double>& vertices_1d_x,
                                      std::vector<double>& vertices_1d_y,
                                      std::vector<double>& vertices_1d_z);
}

#include "chi_meshvector.h"
#include "chi_meshmatrix3x3.h"
#include "chi_meshtensor_rank2_dim3.h"
#include "chi_meshface.h"
#include "chi_mesh_edgeloops.h"
#include "chi_mesh_interface.h"

#include "ChiMesh/SweepUtilities/sweep_namespace.h"



#endif