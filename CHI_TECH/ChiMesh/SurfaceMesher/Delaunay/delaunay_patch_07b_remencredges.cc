#include "delaunay_mesher.h"

//###################################################################
/**Checks for and removes encroached edges.*/
bool chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::
RemoveEncroachedEdges()
{
  bool encroachedEdgeFound = false;

  for (unsigned t=0; t<triangles.size(); t++)
  {
    if (!triangles[t].invalidated)
    {
      //printf("Tri %d\n",t);
      for (int e=0; e<3; e++)
      {
        //printf("Edge %d\n",e);
        if (triangles[t].e_index[e][2] < 0)
        {
          if (CheckEdgeCanBeSplit(triangles[t].e_index[e][0],
                                  triangles[t].e_index[e][1]))
          {
            int i_index = triangles[t].e_index[e][0];
            int f_index = triangles[t].e_index[e][1];
            int o_index;
            if (e == 2)
            {
              o_index = triangles[t].e_index[0][1];
            }
            else
            {
              o_index = triangles[t].e_index[e+1][1];
            }

            if (CheckEdgeEncroached(i_index,f_index,o_index))
            {
              Vertex v0 = Pstar[i_index];
              Vertex v1 = Pstar[f_index];

              Vertex vc = v0*0.5 + v1*0.5;

              printf("Encroached edge %d->%d\n",i_index,f_index);
              if (InsertVertex(vc, t))
              {
                encroachedEdgeFound = true;
                break;
              }

            }
          }
        }//if edge<0
        if (encroachedEdgeFound) break;
      }//for edge e
      if (encroachedEdgeFound) break;
    }

  }//for triangle t

  return encroachedEdgeFound;
}