#ifndef CHI_MESH_MESH_CUTTING_H
#define CHI_MESH_MESH_CUTTING_H

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/Cell/cell_polygon.h"
#include "ChiMesh/Cell/cell_polyhedron.h"

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
    uint64_t cut_point_id=0;

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

  //2D_utils
  std::pair<uint64_t,uint64_t>
    MakeEdgeFromPolygonEdgeIndex(const std::vector<uint64_t>& vertex_ids,
                                 int edge_index);

  void PopulatePolygonFromVertices(
    const MeshContinuum& mesh,
    const std::vector<uint64_t>& vertex_ids,
    chi_mesh::Cell& cell);

  bool CheckPolygonQuality(const MeshContinuum& mesh, const chi_mesh::Cell& cell);

  void SplitConcavePolygonsIntoTriangles(MeshContinuum& mesh,
                                         std::vector<chi_mesh::Cell*>& cell_list);

  void CutPolygon(const std::vector<ECI>& cut_edges,
                  const std::set<uint64_t>& cut_vertices,
                  const Vector3 &plane_point,
                  const Vector3 &plane_normal,
                  MeshContinuum& mesh,
                  chi_mesh::CellPolygon& cell);

  //3D_utils
  bool CheckPolyhedronQuality(const MeshContinuum& mesh,
                              const chi_mesh::Cell& cell);

  void PopulatePolyhedronFromFaces(
    const MeshContinuum& mesh,
    const std::vector<std::vector<uint64_t>>& raw_faces,
    chi_mesh::Cell& cell);

  void CutTetrahedron(const std::vector<ECI>& global_cut_edges,
                      const std::set<uint64_t>& global_cut_vertices,
                      const Vector3 &plane_point,
                      const Vector3 &plane_normal,
                      MeshContinuum& mesh,
                      chi_mesh::CellPolyhedron& cell);

  //plane
  void CutMeshWithPlane(MeshContinuum& mesh,
                        const Vector3& plane_point,
                        const Vector3& plane_normal,
                        double merge_tolerance=1.0e-3,
                        double float_compare=1.0e-10);
}
}

#endif //CHI_MESH_MESH_CUTTING_H