#include "delaunay_mesher.h"

//###################################################################
/**Setup first triangle*/
void chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::CreateFirstTriangle()
{
  //================================================== First hull point
  int hi0=-1; //Indices
  int hi1=-1;
  int hi2=-1;

  Vertex hv0; //Vertices
  Vertex hv1;
  Vertex hv2;

  std::vector<int>::reverse_iterator hull_iter;

  //================================================== First hull point
  hi0 = this->unused_vs.back(); this->unused_vs.pop_back();
  hv0 = Pstar.at(hi0);

  //================================================== Second hull point
  for (hull_iter = this->unused_vs.rbegin();
       hull_iter != this->unused_vs.rend();
       hull_iter++)
  {
    Vertex hullv = Pstar.at(*hull_iter);

    if (fabs(hullv.x-hv0.x)>tolerance)
    {
      hi1 = *hull_iter;
      hv1 = hullv;
      this->unused_vs.erase((hull_iter+1).base());
      break;
    }
  }

  //================================================== Third hull point
  hi2 = this->unused_vs.back(); this->unused_vs.pop_back();
  hv2 = Pstar.at(hi2);

  int a,b,c;
  if (this->Orient2D(hv0,hv1,hv2)>0.0)
  {
    a = hi0;
    b = hi1;
    c = hi2;
  }
  else
  {
    a = hi0;
    b = hi2;
    c = hi1;
  }

  //================================================== Create first triangle
  Tri first_tri;
  first_tri.SetIndices(a,b,c);
  first_tri.geometric_normal = this->patch_normal;

  //================================================== InitializeAlphaElements open edges
  Edge e1,e2,e3;
  e1.v_index[0] = first_tri.v_index[0];
  e1.v_index[1] = first_tri.v_index[1];
  e1.f_index[0] = 0;
  e1.f_index[1] = -1;
  e1.f_index[2] = 0;
  e1.vertices[0] = hv0;
  e1.vertices[1] = hv1;

  e2.v_index[0] = first_tri.v_index[1];
  e2.v_index[1] = first_tri.v_index[2];
  e2.f_index[0] = 0;
  e2.f_index[1] = -1;
  e1.f_index[2] = 1;
  e2.vertices[0] = hv1;
  e2.vertices[1] = hv2;

  e3.v_index[0] = first_tri.v_index[2];
  e3.v_index[1] = first_tri.v_index[0];
  e3.f_index[0] = 0;
  e3.f_index[1] = -1;
  e1.f_index[2] = 2;
  e3.vertices[0] = hv2;
  e3.vertices[1] = hv0;

  this->open_edges.push_back(e1);
  this->open_edges.push_back(e2);
  this->open_edges.push_back(e3);


  //printf("First Tri = %d->%d->%d\n",a,b,c);

  this->triangles.push_back(first_tri);
}