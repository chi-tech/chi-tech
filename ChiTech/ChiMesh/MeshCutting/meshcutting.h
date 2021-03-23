#ifndef CHI_MESH_MESH_CUTTING_H
#define CHI_MESH_MESH_CUTTING_H

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/Cell/cell_polygon.h"

namespace chi_mesh
{
namespace mesh_cutting
{
  typedef std::pair<uint64_t, uint64_t> Edge;

  /**A general structure storing the two vertex-ids
   * to be cut or being cut. Additionally a cut-point
   * id is stored.*/
  struct ECI    //EdgeCutInfo
  {
    Edge vertex_ids;
    uint64_t cut_point_id;

    ECI() = default;

    explicit ECI(Edge in_edge,
                 uint64_t in_cutpoint_id) :
      vertex_ids(std::move(in_edge)),
      cut_point_id(in_cutpoint_id)
    {}

    static
    bool Comparator(const ECI& edge_cut_info,
                    const Edge& ref_edge)
    {
      return ref_edge == edge_cut_info.vertex_ids;
    }
  };

void CutMeshWithPlane(MeshContinuum& mesh,
                      const Vector3& plane_point,
                      const Vector3& plane_normal);
void CutPolygon(const std::vector<ECI>& cut_edges,
                const std::set<uint64_t>& cut_vertices,
                MeshContinuum& mesh,
                chi_mesh::CellPolygon& cell);
}
}

#endif //CHI_MESH_MESH_CUTTING_H