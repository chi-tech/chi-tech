#include "delaunay_mesher.h"
#include <unistd.h>

//###################################################################
/**Attaches unused vertices to a convex hull.*/
void chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::
    AttachUnusedVertex(int v_index)
{
  Vertex vc = this->Pstar.at(v_index); unused_vs.pop_back();
  bool attachment_made = false;
  bool intersect_simplex = false;

  //================================================== Loop over open edges
  std::vector<Edge>::iterator curedge;
  for (curedge = open_edges.begin(); curedge!=open_edges.end(); curedge++)
  {
    Vertex va = curedge->vertices[0];
    Vertex vb = curedge->vertices[1];

    if (Orient2D(va,vb,vc)<(-tolerance))
    {
      bool bad_leg0 = CheckCrossSimplices(va,vc);
      bool bad_leg1 = CheckCrossSimplices(vc,vb);
      bad_leg0 = false;
      bad_leg1 = false;

      if ((!bad_leg0) && (!bad_leg1))
      {
        //======================================= Create new triangle
        int va_index = curedge->v_index[0];
        int vb_index = curedge->v_index[1];

        Tri new_tri;
        new_tri.SetIndices(va_index,v_index,vb_index);
        new_tri.geometric_normal = this->patch_normal;
        new_tri.e_index[2][2] = curedge->f_index[0];
        new_tri.e_index[2][3] = curedge->f_index[2];

        this->triangles.push_back(new_tri);

        //======================================= Update open edges
        Edge new_edge1;
        new_edge1.v_index[0] = va_index;
        new_edge1.v_index[1] = v_index;
        new_edge1.vertices[0]= va;
        new_edge1.vertices[1]= vc;
        new_edge1.f_index[0] = this->triangles.size()-1;
        new_edge1.f_index[2] = 0;

        Edge new_edge2;
        new_edge2.v_index[0] = v_index;
        new_edge2.v_index[1] = vb_index;
        new_edge2.vertices[0]= vc;
        new_edge2.vertices[1]= vb;
        new_edge2.f_index[0] = this->triangles.size()-1;
        new_edge2.f_index[2] = 1;

        open_edges.erase(curedge);
        open_edges.insert(curedge,new_edge2);
        open_edges.insert(curedge,new_edge1);

        printf("VTri added %d->%d->%d\n",va_index,v_index,vb_index);

        attachment_made = true;
        break;
      }
      else
      {
        intersect_simplex = true;
        unused_vs.insert(unused_vs.end()-1,v_index);
        //printf("Vertex kicked back\n");
        //exit(EXIT_FAILURE);
      }
    }


  }

  if ((!attachment_made) && (!intersect_simplex))
  {
    std::cerr << "ERROR: vertex attachment failed.\n";
    exit(EXIT_FAILURE);
  }
}