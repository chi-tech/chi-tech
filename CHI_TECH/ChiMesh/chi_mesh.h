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
  struct Vector;
  typedef Vector Normal;
  typedef Vector Vertex;
  typedef Vector Node;

  struct Matrix3x3;

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
  class CellSlab;
  class CellTriangle;
  class CellPolygon;
  class CellPolyhedron;

  //=================================== Field function interpolation
  class FieldFunctionInterpolation;
  class FieldFunctionInterpolationSlice;
  class FieldFunctionInterpolationLine;
  class FieldFunctionInterpolationVolume;

  //=================================== Meshes
  class LineMesh;
  class SurfaceMesh;
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
  class SurfaceMesherPredefined;
  class SurfaceMesherDelaunay;
  class SurfaceMesherTriangle;
  //struct DelaunayMeshContext;

  struct Interface;

  //==================================== Volume meshers
  class VolumeMesher;
  class VolumeMesherLinemesh1D;
  class VolumeMesherPredefined2D;
  class VolumeMesherExtruder;




  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ROUTINES

  Boundary*              AssignSurfaceToBoundary(chi_mesh::SurfaceMesh* surface);
  MeshHandler*           GetCurrentHandler();
  EdgeLoopCollection*    SplitEdgeLoopByAngle(EdgeLoop* input,double angle=1);

  //=================================== Domain decompositions
  double ComputeLBF(std::vector<Vector>& points,
                    std::vector<double>& x_cuts,
                    std::vector<double>& y_cuts);
  void   DecomposeSurfaceMeshPxPy(SurfaceMesh* smesh, int Px, int Py);
  void   SurfaceMeshImprintLine(SurfaceMesh* smesh,
                                Vertex* line_point_0,
                                Vertex* line_point_1,
                                double tolerance);




  //=================================== Raytracing
  void RayTrace(chi_mesh::MeshContinuum* grid,
                Cell* cell,
                const Vector& pos_i,
                const Vector& omega_i,
                double& d_to_surface,
                Vector& pos_f,
                int* aux_info = nullptr);
  bool
  CheckPlaneLineIntersect(chi_mesh::Normal plane_normal,
                          chi_mesh::Vector plane_point,
                          chi_mesh::Vector line_point_0,
                          chi_mesh::Vector line_point_1,
                          chi_mesh::Vector& intersection_point,
                          std::pair<double,double>& weights);
}

#include"chi_meshvector.h"
#include "chi_meshmatrix3x3.h"
#include"chi_meshface.h"
#include"chi_mesh_edgeloops.h"
#include"chi_mesh_interface.h"

#include "SweepUtilities/chi_sweep.h"



#endif