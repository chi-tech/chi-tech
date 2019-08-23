#include "delaunay_mesher.h"

//###################################################################
/**Checks for and removes bad quality triangles.*/
bool chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::
RemoveBadQualityTriangles()
{
  bool bad_triangles_found = false;

  for (unsigned t=0; t<triangles.size(); t++)
  {
    if (!triangles[t].invalidated)
    {
      int v0_index = triangles[t].v_index[0];
      int v1_index = triangles[t].v_index[1];
      int v2_index = triangles[t].v_index[2];

      Vertex v0 = Pstar[v0_index];
      Vertex v1 = Pstar[v1_index];
      Vertex v2 = Pstar[v2_index];

      //========================================= Find smallest edge
      Vector v01 = v1 - v0;
      Vector v12 = v2 - v1;
      Vector v20 = v0 - v2;

      double d01 = v01.Norm();
      double d12 = v12.Norm();
      double d20 = v20.Norm();

      double avg_edge_length = (d01+d12+d20)/3.0;

      double smallest_edge = d01;

      if (d12<smallest_edge) smallest_edge = d12;
      if (d20<smallest_edge) smallest_edge = d20;

      //========================================= Find circum-circle
      Vertex circum_center;
      double radius;

      FindCircumCircle(v0,v1,v2,circum_center,radius);



      //========================================= Check quality
      if ( ((radius/smallest_edge)>quality_factor1) ||
           (avg_edge_length>(2.0*absoluteMinumumSize)) )
      {
        //printf("Quality evaluated for triangle %d = %f ",t,radius/smallest_edge);

        bool has_encroached_edge=false;
        bool causes_encroached_edge=false;
        for (int e=0; e<3; e++)
        {
          if (triangles[t].e_index[e][2]<0)
          {
            //================================== Check if it has an encroached edge
            if (e==2)
            {
              has_encroached_edge =
                CheckEdgeEncroached(triangles[t].e_index[e][0],
                                    triangles[t].e_index[e][1],
                                    triangles[t].e_index[0][1]);
            } else
            {
              has_encroached_edge =
                CheckEdgeEncroached(triangles[t].e_index[e][0],
                                    triangles[t].e_index[e][1],
                                    triangles[t].e_index[e+1][1]);
            }

            //================================== Check if it will cause an encroached edge
            if (triangles[t].e_index[e][2]<0)
            {
              causes_encroached_edge =
                CheckEdgeEncroached(triangles[t].e_index[e][0],
                                    triangles[t].e_index[e][1],
                                    circum_center);
            }
          }





           if (causes_encroached_edge || has_encroached_edge) break;
        }

        if ((!has_encroached_edge) && (!causes_encroached_edge))
        {
          //printf("Removing it center=%.3f,%.3f\n",circum_center.x, circum_center.y);
          if (InsertVertex(circum_center, t))
          {
            bad_triangles_found = true;
            break;
          }

        }
        else
        {
          //printf("Encroached and therefore will not be removed\n");
        }

      }

    }
  }//for t

  return bad_triangles_found;
}