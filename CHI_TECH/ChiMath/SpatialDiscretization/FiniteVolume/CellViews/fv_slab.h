#ifndef _fv_slab_h
#define _fv_slab_h

#include "fv_cellbase.h"
#include <CHI_MESH/CHI_CELL/cell_slab.h>
#include <CHI_MESH/CHI_MESHCONTINUUM/chi_meshcontinuum.h>


//######################################################### Class def
/**Finite Volume implementation for a slab.*/
class SlabFVView : public CellFVView
{
private:
  chi_mesh::MeshContinuum* grid;
  int v0i;
  int v1i;

public:
  double volume; ///< Actually length times unity dx*dy

  SlabFVView(chi_mesh::CellSlab *slab_cell,
             chi_mesh::MeshContinuum *vol_continuum) :
             CellFVView(2)
  {
    grid = vol_continuum;
    v0i = slab_cell->v_indices[0];
    v1i = slab_cell->v_indices[1];
    chi_mesh::Vertex v0 = *grid->nodes[v0i];
    chi_mesh::Vertex v1 = *grid->nodes[v1i];

    chi_mesh::Vector v01 = v1-v0;
    volume = v01.Norm();
  }
};

#endif