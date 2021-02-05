#ifndef _chi_mesh_cell_polygon_h
#define _chi_mesh_cell_polygon_h

#include "cell.h"
#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>

namespace chi_mesh
{

//######################################################### Class def
/** Polygon cell definition.*/
class CellPolygon : public Cell
{
private:
  std::vector<chi_mesh::Vector3> segment_normals;
  bool segment_normals_developed = false;
public:
  CellPolygon() : Cell(CellType::POLYGON) {}

private:
  void DevelopSegmentNormals(const chi_mesh::MeshContinuum& grid)
  {
    segment_normals.reserve(faces.size());
    for (auto& face : faces) //edges
    {
      chi_mesh::Vertex &v0 = *grid.vertices[face.vertex_ids[0]];
      const chi_mesh::Vertex &vc = centroid;

      chi_mesh::Vector3 khat(0.0, 0.0, 1.0);
      auto vc0 = vc - v0;
      auto n0 = ((vc0) / vc0.Norm()).Cross(khat);
      n0 = n0 / n0.Norm();

      segment_normals.push_back(n0);
    }
    segment_normals_developed = true;
  }
public:
  std::vector<chi_mesh::Vector3>& GetSegmentNormals(const chi_mesh::MeshContinuum& grid)
  {
    if (!segment_normals_developed)
      DevelopSegmentNormals(grid);

    return segment_normals;
  }
};

}

#endif