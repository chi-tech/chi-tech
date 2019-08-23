#include "cell_triangle.h"
#include <iostream>

//###################################################################
/**Checks if an unconnected edge of a cell has a boundary condition
 * applied to it.

 \return Returns the amount of unconnected edges.
 * */
bool chi_mesh::CellTriangle::CheckBoundary2D()
{
  for (int e=0;e<3;e++)
  {
    if ((this->e_index[e][2]<0) && (this->e_index[e][3]<0))
    {
      return false;
    }
  }

  return true;
}