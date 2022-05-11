#ifndef CHI_MESHHANDLER_H
#define CHI_MESHHANDLER_H

#include <vector>

#include"../chi_mesh.h"



/**Object for containing all mesh related variables.*/
class chi_mesh::MeshHandler
{
public:
  std::shared_ptr<chi_mesh::SurfaceMesher> surface_mesher = nullptr;
  std::shared_ptr<chi_mesh::VolumeMesher>  volume_mesher = nullptr;


public:
  chi_mesh::MeshContinuumPtr& GetGrid() const;
  MeshHandler() = default;
  MeshHandler(const MeshHandler&) = delete;
  MeshHandler& operator=(const MeshHandler&) = delete;
};

#endif//CHI_MESHHANDLER_H