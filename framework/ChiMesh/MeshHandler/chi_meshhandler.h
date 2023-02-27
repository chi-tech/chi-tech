#ifndef CHI_MESHHANDLER_H
#define CHI_MESHHANDLER_H

#include <utility>
#include <vector>

#include"../chi_mesh.h"



/**Object for containing all mesh related variables.*/
class chi_mesh::MeshHandler
{
protected:
  std::shared_ptr<chi_mesh::SurfaceMesher> surface_mesher_ = nullptr;
  std::shared_ptr<chi_mesh::VolumeMesher>  volume_mesher_ = nullptr;

public:
  chi_mesh::MeshContinuumPtr& GetGrid() const;
  chi_mesh::SurfaceMesher& GetSurfaceMesher();
  chi_mesh::VolumeMesher& GetVolumeMesher();
  const chi_mesh::SurfaceMesher& GetSurfaceMesher() const;
  const chi_mesh::VolumeMesher& GetVolumeMesher() const;
  void SetSurfaceMesher(std::shared_ptr<chi_mesh::SurfaceMesher> surface_mesher)
  {surface_mesher_ = std::move(surface_mesher);}
  void SetVolumeMesher(std::shared_ptr<chi_mesh::VolumeMesher> volume_mesher)
  {volume_mesher_ = std::move(volume_mesher);}
  MeshHandler() = default;
  MeshHandler(const MeshHandler&) = delete;
  MeshHandler& operator=(const MeshHandler&) = delete;
};

#endif//CHI_MESHHANDLER_H