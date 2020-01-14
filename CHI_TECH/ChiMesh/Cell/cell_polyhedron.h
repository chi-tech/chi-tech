#ifndef _chi_mesh_cell_polyhedron_h
#define _chi_mesh_cell_polyhedron_h

#include "cell.h"
#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>

namespace chi_mesh
{


//######################################################### Class def
/** Polyhedron cell definition.*/
class CellPolyhedron : public Cell
{
private:
  std::vector<std::vector<std::vector<int>>> face_edges;
  bool edges_developed = false;

private:
  std::vector<std::vector<chi_mesh::Vector>> face_segment_normals;
  bool segment_normals_developed = false;

public:
  CellPolyhedron() : Cell(CellType::POLYHEDRON) {}

private:
  void DevelopEdges()
  {
    size_t num_faces = faces.size();
    for (size_t f=0; f<num_faces; f++)
    {
      CellFace& face = faces[f];

      typedef std::vector<int> VecInt;
      typedef std::vector<VecInt> VecVecInt;

      size_t num_face_verts = face.vertex_ids.size();
      VecVecInt edges(num_face_verts,VecInt(2,-1));
      for (size_t fvi=0; fvi<num_face_verts; fvi++)
      {
        edges[fvi][0] = face.vertex_ids[fvi];
        if (fvi < (num_face_verts-1))
          edges[fvi][1] = face.vertex_ids[fvi+1];
        else
          edges[fvi][1] = face.vertex_ids[0];
      }
      face_edges.push_back(edges);
    }//for f
    edges_developed = true;
  }

  void DevelopSegmentNormals(const chi_mesh::MeshContinuum* grid)
  {
    auto& vcc = centroid;

    face_segment_normals.resize(faces.size());
    int f=-1;
    for (auto& face : faces)
    {
      f++;
      auto& vfc = face.centroid;
      face_segment_normals[f].reserve(face.vertex_ids.size());
      for (auto vi : face.vertex_ids)
      {
        auto& vert = *grid->nodes[vi];

        auto v01 = vfc - vert;
        auto v12 = vcc - vfc;

        auto n = v01.Cross(v12);
        n = n/n.Norm();

        face_segment_normals[f].push_back(n);
      }
    }
    segment_normals_developed = true;
  }

public:
  std::vector<std::vector<int>>& GetFaceEdges(size_t f)
  {
    if (!edges_developed)
      DevelopEdges();
    return face_edges[f];
  }

  std::vector<std::vector<chi_mesh::Vector>>&
    GetSegmentNormals(const chi_mesh::MeshContinuum* grid)
  {
    if (!segment_normals_developed)
      DevelopSegmentNormals(grid);

    return face_segment_normals;
  }

};

}

#endif