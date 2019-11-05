#ifndef _chi_mesh_cell_slab_h
#define _chi_mesh_cell_slab_h

#include "cell_newbase.h"

namespace chi_mesh
{

//######################################################### Class def
/** Slab cell definition.*/
class CellSlabV2 : public CellBase
{
public:
  CellSlabV2() : CellBase(CellType::SLABV2) {}


};

}

#endif