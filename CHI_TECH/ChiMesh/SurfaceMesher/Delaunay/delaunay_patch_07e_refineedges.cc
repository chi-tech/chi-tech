#include "delaunay_mesher.h"

//###################################################################
/**Checks for and removes bad quality triangles.*/
bool chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::
RefineEdges()
{
  bool edge_split=false;
  for (unsigned t=0; t<triangles.size(); t++)
  {
    if (!triangles[t].invalidated)
    {
      for (int e=0; e<3; e++)
      {
        if (triangles[t].e_index[e][2]>=0)
        {
          int v0_i = triangles[t].e_index[e][0];
          int v1_i = triangles[t].e_index[e][1];

          Vertex v0 = Pstar[v0_i];
          Vertex v1 = Pstar[v1_i];

          Vector3 v01 = v1 - v0;
          double v01_ell = v01.Norm();

          if (v01_ell>(2.0*absoluteMinumumSize))
          {
            Vertex edge_center = v0 + v01*0.5;

            if (InsertVertex(edge_center, t))
            {
              edge_split = true;
              break; //from for e
            }

          }
        }
      }
    }
    if (edge_split) break; //from for t
  }//for t

  return edge_split;
}