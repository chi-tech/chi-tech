#ifndef _cell_slab_h
#define _cell_slab_h

#include "cell.h"

//################################################################### Class def
/**Object to handle one dimensional slab cells.
 *
 * edges\n
 * [0] Neighbor cell index to the left.
 * [1] Neighbor cell index to the right.*/
class chi_mesh::CellSlab : public chi_mesh::Cell
{
public:
  int     v_indices[2];
  int     edges[2];
  Vector  face_normals[2];
};

#endif
