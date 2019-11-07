#ifndef _chi_mesh_cell_polyhedron_h
#define _chi_mesh_cell_polyhedron_h

#include "cell.h"

namespace chi_mesh
{


//######################################################### Class def
/** Polyhedron cell definition.*/
class CellPolyhedronV2 : public Cell
{
public:
  CellPolyhedronV2() : Cell(CellType::POLYHEDRONV2) {}

  std::vector<std::vector<int>> GetFaceEdges(size_t f)
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

    return edges;
  }
};

}

#endif