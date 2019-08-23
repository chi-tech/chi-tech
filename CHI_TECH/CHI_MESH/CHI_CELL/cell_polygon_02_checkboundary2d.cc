#include "cell_polygon.h"
#include <iostream>

//###################################################################
/**Checks if an unconnected edge of a cell has a boundary condition
 * applied to it.

 * */
bool chi_mesh::CellPolygon::CheckBoundary2D()
{
  for (int e=0;e<edges.size();e++)
  {
    if ((this->edges[e][2]<0) && (this->edges[e][3]<0))
    {
      return false;
    }
  }

  return true;
}