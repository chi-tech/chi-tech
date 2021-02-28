#ifndef CHI_MESH_CELL_SLAB_H
#define CHI_MESH_CELL_SLAB_H

#include "cell.h"

namespace chi_mesh
{

//######################################################### Class def
/** Slab cell definition.*/
class CellSlab : public Cell
{
public:
  CellSlab() : Cell(CellType::SLAB) {}
};

}

#endif