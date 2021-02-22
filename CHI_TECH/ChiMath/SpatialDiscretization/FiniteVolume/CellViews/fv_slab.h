#ifndef SLAB_FV_VALUES_H
#define SLAB_FV_VALUES_H

#include "fv_cellbase.h"
#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>


//######################################################### Class def
/**Finite Volume implementation for a slab.*/
class SlabFVValues : public CellFVValues
{
private:
  chi_mesh::MeshContinuumPtr grid;
  int v0i;
  int v1i;

public:

  SlabFVValues(chi_mesh::CellSlab *slab_cell,
               chi_mesh::MeshContinuumPtr vol_continuum) :
    CellFVValues(2)
  {
    grid = vol_continuum;
    v0i = slab_cell->vertex_ids[0];
    v1i = slab_cell->vertex_ids[1];
    chi_mesh::Vertex v0 = *grid->vertices[v0i];
    chi_mesh::Vertex v1 = *grid->vertices[v1i];

    chi_mesh::Vector3 v01 = v1 - v0;
    volume = v01.Norm();
    face_area.push_back(1.0);
    face_area.push_back(1.0);
  }
};

#endif