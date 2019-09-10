#ifndef _cell_polygon_h
#define _cell_polygon_h

#include "cell.h"
#include"../Region/chi_region.h"

#define EDGE_NEIGHBOR 2
//######################################################### Class def
/**Object to handle generic polygon cells.
 *
 * edges\n
 An array of 4 integers.\n
 [0] = Vertex index of edge start.\n
 [1] = Vertex index of edge end.\n
 [2] = Index of the cell adjoining this edge (not the current cell).
       -1 if not connected to anything,-1*boundary_index if connected
       to a boundary.\n
 [3] = Edge number of adjoining face. -1 if not connected
       to anything. 0 if a boundary.\n
   \n*/
class chi_mesh::CellPolygon : public chi_mesh::Cell
{
public:
  std::vector<int>  v_indices;
  std::vector<int*> edges; ///< Stores arrays of edge indices
  std::vector<chi_mesh::Vector> edgenormals;

  //01
  void FindBoundary2D(chi_mesh::Region* region);
  //02
  bool CheckBoundary2D();
};

#endif