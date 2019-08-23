#include "delaunay_mesher.h"

//###################################################################
/***/
bool chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::ConvexifyHull()
{
  bool modification_made =false;
  //printf("Convexifying hull\n");
  //================================================== Loop over open edges
  std::vector<Edge>::iterator curedge;
  for (curedge = open_edges.begin(); curedge!=(open_edges.end()-1); curedge++)
  {
    Vertex va = curedge->vertices[0];
    Vertex vb = curedge->vertices[1];
    Vertex vc = (curedge+1)->vertices[1];

    if (Orient2D(va,vb,vc)<(-tolerance))
    {
      bool bad_leg0 = CheckCrossSimplices(va,vc);
      bad_leg0 =false;
      if (!bad_leg0)
      {
        //======================================= Create new triangle
        int va_index = curedge->v_index[0];
        int vb_index = curedge->v_index[1];
        int vc_index = (curedge+1)->v_index[1];

        Tri new_tri;
        new_tri.SetIndices(va_index,vc_index,vb_index);
        new_tri.geometric_normal = this->patch_normal;

        new_tri.e_index[1][2] = (curedge+1)->f_index[0];
        new_tri.e_index[1][3] = (curedge+1)->f_index[2];

        new_tri.e_index[2][2] = curedge->f_index[0];
        new_tri.e_index[2][3] = curedge->f_index[2];

        this->triangles.push_back(new_tri);

        //======================================= Update open edges
        Edge new_edge1;
        new_edge1.v_index[0] = va_index;
        new_edge1.v_index[1] = vc_index;
        new_edge1.vertices[0]= va;
        new_edge1.vertices[1]= vc;
        new_edge1.f_index[0] = this->triangles.size()-1;
        new_edge1.f_index[2] = 0;


        open_edges.erase(curedge);
        open_edges.erase(curedge);
        open_edges.insert(curedge,new_edge1);



        printf("CTri = %d->%d->%d\n",va_index,vc_index,vb_index);
        modification_made=true;
        break;
      }
      else
      {
        break;
      }
    }
  }
  return modification_made;
}