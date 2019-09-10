#include"delaunay_mesher.h"

//###################################################################
/**Restores the Delaunay conditioning of the mesh using edge-flipping.
 **/
void chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::RestoreDelaunay()
{
  //================================================= Find initial non-locally
  //                                                  delaunay edge
  std::vector<Edge>   non_loc_del_edges;
  FindNonLocDelEdge(non_loc_del_edges);

  //DumpToScilab("ZScilabMesh.sce");

  //================================================= Iterate to remove
  //                                                  non-locally delaunay edges
  int iter=0;
  while ((non_loc_del_edges.size()>0))
  {
    iter++;

    EdgeFlip(non_loc_del_edges.back());
    non_loc_del_edges.pop_back();


    FindNonLocDelEdge(non_loc_del_edges);
    //if (iter>2) break;

  }



}

