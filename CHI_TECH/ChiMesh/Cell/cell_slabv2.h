#ifndef _chi_mesh_cell_slab_h
#define _chi_mesh_cell_slab_h

#include "cell.h"

namespace chi_mesh
{

//######################################################### Class def
/** Slab cell definition.*/
class CellSlabV2 : public Cell
{
public:
  CellSlabV2() : Cell(CellType::SLABV2) {}


};

}

#endif