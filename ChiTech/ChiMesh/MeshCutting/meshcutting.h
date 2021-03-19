#ifndef CHI_MESH_MESH_CUTTING_H
#define CHI_MESH_MESH_CUTTING_H

#include "ChiMesh/chi_mesh.h"

namespace chi_mesh
{
namespace mesh_cutting
{
void CutMeshWithPlane(MeshContinuum& mesh,
                      const Vector3& plane_point,
                      const Vector3& plane_normal);
}
}

#endif //CHI_MESH_MESH_CUTTING_H