#ifndef SLAB_FV_VALUES_H
#define SLAB_FV_VALUES_H

#include "fv_cellbase.h"

//######################################################### Class def
namespace chi_math
{
  /**Finite Volume implementation for a slab.*/
  class SlabFVValues : public chi_math::CellFVValues
  {
  private:
    uint64_t v0i;
    uint64_t v1i;

  public:

    SlabFVValues(const chi_mesh::Cell &slab_cell,
                 const chi_mesh::MeshContinuumConstPtr &grid);
  };
}

#endif