#ifndef _chi_mesh_cell_polygon_h
#define _chi_mesh_cell_polygon_h

#include "cell_newbase.h"

namespace chi_mesh
{

//######################################################### Class def
/** Polygon cell definition.*/
class CellPolygonV2 : public CellBase
{
public:
  CellPolygonV2() : CellBase(CellType::POLYGONV2) {}

  void FindBoundary2D(chi_mesh::Region* region);
  //02
  bool CheckBoundary2D();
};

}

#endif