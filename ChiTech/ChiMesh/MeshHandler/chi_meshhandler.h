#ifndef CHI_MESHHANDLER_H
#define CHI_MESHHANDLER_H

#include <vector>

#include"../chi_mesh.h"


/**Object for containing all mesh related variables.*/
class chi_mesh::MeshHandler
{
public:
  std::vector<chi_mesh::SurfaceMesh*>                surface_mesh_stack;
  std::vector<chi_mesh::LogicalVolume*>              logicvolume_stack;
  std::vector<chi_mesh::FieldFunctionInterpolation*> ffinterpolation_stack;
  std::vector<chi_mesh::UnpartitionedMesh*>          unpartitionedmesh_stack;

  chi_mesh::SurfaceMesher* surface_mesher = nullptr;
  chi_mesh::VolumeMesher*  volume_mesher = nullptr;


public:
  chi_mesh::MeshContinuumPtr& GetGrid() const;
  MeshHandler() = default;
  MeshHandler(const MeshHandler&) = delete;
  MeshHandler& operator=(const MeshHandler&) = delete;
};

#endif//CHI_MESHHANDLER_H