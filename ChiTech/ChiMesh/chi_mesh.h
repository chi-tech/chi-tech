#ifndef CHI_MESH_H
#define CHI_MESH_H

#include<vector>
#include<iostream>
#include<memory>


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

  struct Matrix3x3;
  struct TensorRank2Dim3;

  struct Face;
  struct Edge;
  struct PolyFace;


  struct SPDS;

  //=================================== Cells
  class Cell;

  //=================================== Field function interpolation
  class FieldFunctionInterpolation;
  class FieldFunctionInterpolationSlice;
  class FieldFunctionInterpolationLine;
  class FieldFunctionInterpolationVolume;

  //=================================== Meshes
  class SurfaceMesh;
  class UnpartitionedMesh;
  class MeshContinuum;
  typedef std::shared_ptr<MeshContinuum> MeshContinuumPtr;

  //=================================== Logical Volumes
  class LogicalVolume;
  class SphereLogicalVolume;
  class RPPLogicalVolume;
  class RCCLogicalVolume;
  class SurfaceMeshLogicalVolume;
  class BooleanLogicalVolume;

  //=================================== Mesh handler
  class MeshHandler;

  //=================================== Surface Meshers
  class SurfaceMesher;
  class SurfaceMesherPassthrough;
  class SurfaceMesherPredefined;

  //==================================== Volume meshers
  class VolumeMesher;
  class VolumeMesherExtruder;
  class VolumeMesherPredefinedUnpartitioned;




  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ROUTINES
  MeshHandler&           GetCurrentHandler();
  size_t                 PushNewHandlerAndGetIndex();

  //=================================== Domain decompositions
  double ComputeLBF(std::vector<Vector3>& points,
                    std::vector<double>& x_cuts,
                    std::vector<double>& y_cuts);
  void   DecomposeSurfaceMeshPxPy(const SurfaceMesh& smesh, int Px, int Py);

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

#include "ChiMesh/SweepUtilities/sweep_namespace.h"



#endif