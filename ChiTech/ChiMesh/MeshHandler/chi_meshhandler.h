#ifndef CHI_MESHHANDLER_H
#define CHI_MESHHANDLER_H

#include <vector>

#include"../chi_mesh.h"
#include "../LineMesh/chi_linemesh.h"


/**Object for containing all mesh related variables.*/
class chi_mesh::MeshHandler
{
public:
  typedef std::vector<chi_mesh::SurfaceMesh*> SurfaceMeshCollection;
  std::vector<SurfaceMeshCollection*>                surface_mesh_collections;
  std::vector<chi_mesh::SurfaceMesh*>                surface_mesh_stack;
  std::vector<chi_mesh::LineMesh*>                   linemesh_stack;
  std::vector<chi_mesh::LogicalVolume*>              logicvolume_stack;
  std::vector<chi_mesh::FieldFunctionInterpolation*> ffinterpolation_stack;
  std::vector<chi_mesh::UnpartitionedMesh*>          unpartitionedmesh_stack;

  std::vector<chi_mesh::EdgeLoopCollection*>         edge_loop_collections;

  chi_mesh::SurfaceMesher* surface_mesher = nullptr;
  chi_mesh::VolumeMesher*  volume_mesher = nullptr;


public:
  chi_mesh::MeshContinuumPtr& GetGrid() const;
  MeshHandler() = default;
  MeshHandler(const MeshHandler&) = delete;
  MeshHandler& operator=(const MeshHandler&) = delete;
};

#endif//CHI_MESHHANDLER_H