#ifndef _chi_mesh_cell_polyhedron_h
#define _chi_mesh_cell_polyhedron_h

#include "cell.h"

namespace chi_mesh
{


//######################################################### Class def
/** Polyhedron cell definition.*/
class CellPolyhedronV2 : public Cell
{
private:
  std::vector<std::vector<std::vector<int>>> face_edges;
  bool edges_developed  = false;

public:
  CellPolyhedronV2() : Cell(CellType::POLYHEDRONV2) {}

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

public:
  std::vector<std::vector<int>>& GetFaceEdges(size_t f)
  {
    if (!edges_developed)
      DevelopEdges();
    return face_edges[f];
  }
};

}

#endif