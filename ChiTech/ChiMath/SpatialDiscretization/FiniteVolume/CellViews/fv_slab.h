#ifndef SLAB_FV_VALUES_H
#define SLAB_FV_VALUES_H

#include "fv_cellbase.h"
#include "ChiMesh/Cell/cell.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"


//######################################################### Class def
/**Finite Volume implementation for a slab.*/
class SlabFVValues : public CellFVValues
{
private:
  uint64_t v0i;
  uint64_t v1i;

public:

  SlabFVValues(const chi_mesh::Cell& slab_cell,
               const chi_mesh::MeshContinuum& grid) :
    CellFVValues(2)
  {
    v0i = slab_cell.vertex_ids[0];
    v1i = slab_cell.vertex_ids[1];
    const auto& v0 = grid.vertices[v0i];
    const auto& v1 = grid.vertices[v1i];

    chi_mesh::Vector3 v01 = v1 - v0;
    volume = v01.Norm();
    face_area.push_back(1.0);
    face_area.push_back(1.0);
  }
};

#endif