#ifndef _chi_graph_h
#define _chi_graph_h

#include "../ChiMesh/chi_mesh.h"


/**Cell pointer struct.*/
struct GraphCellInfo
{
  chi_mesh::Cell* cell_ptr;

  GraphCellInfo()
  {
    cell_ptr = nullptr;
  }

  GraphCellInfo& operator=(const GraphCellInfo& other)
  {
    cell_ptr = other.cell_ptr;
    return *this;
  }
};



namespace chi_graph
{
  struct GraphVertex;
  class DirectedGraph;
}


#endif