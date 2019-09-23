#include "cell_polygon.h"
#include <iostream>
#include"../LineMesh/chi_linemesh.h"
#include "../Boundary/chi_boundary.h"

//###################################################################
/**This function operates on cells with an open edge and finds
 * the boundary associated with such an edge. For a 2D operation
 * it only tries to match to the line-mesh associated with the
 * boundary. If a boundary was found then the edge index is set to >=0,
 * i.e. e_index[e][3]=0*/
void chi_mesh::CellPolygon::FindBoundary2D(chi_mesh::Region* region)
{
  //================================================== Check if has boundary
  bool has_boundary=false;
  for (int e=0;e<edges.size();e++)
  {
    if (this->edges[e][2]<0)
    {
      has_boundary = true;
      break;
    }
  }
  if (!has_boundary) { return; }

  //================================================== Use init mesh continuum
  chi_mesh::MeshContinuum* ref_cont = region->volume_mesh_continua.back();

  //================================================== Loop over edges
  for (int e=0;e<edges.size();e++)
  {
    if (this->edges[e][2]<0)
    {
      //================================================ Loop over boundaries
      std::vector<chi_mesh::Boundary*>::iterator bndry;
      for (bndry = region->boundaries.begin();
           bndry != region->boundaries.end();
           bndry++)
      {
        chi_mesh::Vertex v[2];

        try{
          v[0] = *ref_cont->nodes.at(this->edges[e][0]);
          v[1] = *ref_cont->nodes.at(this->edges[e][1]);
        }
        catch (const std::out_of_range& o)
        {
          std::cerr << "Program exception in FindBoundary2D!\n";
          exit(EXIT_FAILURE);
        }
        if ((*bndry)->initial_mesh_continuum.line_mesh!= nullptr)
        {
          //======================================= Get vertices
          chi_mesh::Vertex w[2];

          try{
            w[0] = (*bndry)->initial_mesh_continuum.line_mesh->vertices.front();
            w[1] = (*bndry)->initial_mesh_continuum.line_mesh->vertices.back();
          }
          catch (const std::out_of_range& o)
          {
            std::cerr << "Program exception in FindBoundary2D!\n";
            exit(EXIT_FAILURE);
          }

          //======================================= Check slope
          bool slopes_equal = false;
          chi_mesh::Vector v01 = v[1]-v[0];
          chi_mesh::Vector w01 = w[1]-w[0];

          v01 = v01/v01.Norm();
          w01 = w01/w01.Norm();

          //printf("%.4f %.4f %.4f  ",v01.x,v01.y,v01.z);
          //printf("%.4f %.4f %.4f\n",w01.x,w01.y,w01.z);

          if ((1.0-fabs(v01.Dot(w01)))<0.0000001)
          {
            slopes_equal = true;
          }

          //======================================= Check colinear
          bool colinear = false;
          if (slopes_equal)
          {
            chi_mesh::Vector vw0 = w[0]-v[0];

            if (vw0.Norm()<0.0000001)
            {
              colinear = true;
            }
            else
            {
              vw0 = vw0/vw0.Norm();

              double dot = 1.0-fabs(vw0.Dot(v01));
              //printf("     %.4f %.4f %.4f    %.9f\n",vw0.x,vw0.y,vw0.z,dot);

              if (dot<0.0000001)
              {
                colinear = true;
              }
            }

          }

          //======================================= If both conditions met
          if (slopes_equal && colinear)
          {
            this->edges[e][2] =
              -1*std::distance(region->boundaries.begin(),bndry);
            this->edges[e][3] = 0;
            //printf("Boundary assigned %d\n",this->e_index[e][2]);
          }


        }//if linemesh
      }//for bndry




    }

  }

}
