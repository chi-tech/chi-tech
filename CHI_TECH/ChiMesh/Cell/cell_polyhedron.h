#ifndef _cell_polyhedron_h
#define _cell_polyhedron_h

#include "cell.h"
#include "../Region/chi_region.h"

//################################################################### Class def
/**Object to handle generic polyhedron cells.*/
class chi_mesh::CellPolyhedron : public chi_mesh::Cell
{
public:
  std::vector<int>       v_indices;
  std::vector<PolyFace*> faces;

public:
  CellPolyhedron() : Cell(CellType::POLYHEDRON) {}
};


#endif