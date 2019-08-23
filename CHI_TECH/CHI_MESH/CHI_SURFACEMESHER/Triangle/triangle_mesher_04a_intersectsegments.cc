#include "triangle_mesher.h"
#include "../../CHI_MESHCONTINUUM/chi_meshcontinuum.h"
#include "../../CHI_BOUNDARY/chi_boundary.h"

//###################################################################
/** Computes the intersection of a cut-line and the segments in
 * the referenced surface.*/
void chi_mesh::SurfaceMesherTriangle::
     IntersectSegments(chi_mesh::Line cut_line,
                       chi_mesh::SurfaceMesh *surface_mesh)
{
  //================================================== Precalculate cut_line
  Vertex c0 = cut_line.vertices[0];
  Vertex c1 = cut_line.vertices[1];
  Vertex nc = c1 - c0;

  //================================================== Initialize vertex list
  std::vector<int> inserted_verts;


  //================================================== Loop until no
  //                                                   intersections found
  bool found_intersection = true;
  while (found_intersection)
  {
    found_intersection = false;

    //=========================================== Loop over segments
    for (unsigned ell=0; ell<surface_mesh->lines.size(); ell++)
    {
      //==================================== Calculate segment line
      int    s0_index = surface_mesh->lines[ell].v_index[0];
      int    s1_index = surface_mesh->lines[ell].v_index[1];
      Vertex       s0 = surface_mesh->lines[ell].vertices[0];
      Vertex       s1 = surface_mesh->lines[ell].vertices[1];
      Vertex       ns = s1 - s0;

      double a0,a1,a2,a3,tc,ts;

      //==================================== IF VERTICAL CUT LINE
      if ( fabs(nc.x) < 1.0e-5 ) //Vertical check
      {
        a0 = (s0.y-c0.y)/nc.y;
        a1 = ns.y/nc.y;

        a2 = s0.x - c0.x - a0*nc.x;
        a3 = a1*nc.x - ns.x;

        if ( fabs(a3) > 1.0e-5 ) //parallel check
        {
          ts = a2/a3;
          tc = a0 + a1*ts;

          if ( (ts < (1.0-cut_tol)) && (ts > (0.0+cut_tol)) )
          {
            //printf("Cutting segment (%f,%f)->(%f,%f)\n", s0.x,s0.y,s1.x,s1.y);
            //========================= Erase intersected segment
            surface_mesh->lines.erase(surface_mesh->lines.begin()+ell);

            //========================= Build new lines
            Vertex new_vert = s0 + ns*ts;
            surface_mesh->vertices.push_back(new_vert);
            int new_vert_index = surface_mesh->vertices.size()-1;
            inserted_verts.push_back(new_vert_index);

            //printf("Inserted vertex %d %f,%f\n",new_vert_index,new_vert.x,new_vert.y);

            Line new_line_A, new_line_B;

            new_line_A.vertices[0] = s0;
            new_line_A.vertices[1] = new_vert;

            new_line_A.v_index[0] = s0_index;
            new_line_A.v_index[1] = new_vert_index;

            new_line_B.vertices[0] = new_vert;
            new_line_B.vertices[1] = s1;

            new_line_B.v_index[0] = new_vert_index;
            new_line_B.v_index[1] = s1_index;

            surface_mesh->lines.push_back(new_line_A);
            surface_mesh->lines.push_back(new_line_B);

            found_intersection = true;
            break;

          } //if in tol
        } //if not parallel
      } //if vertical line

      //==================================== IF HORIZONTAL CUT LINE
      if ( fabs(nc.y) < 1.0e-5 ) //Horizontal check
      {
        a0 = (s0.x-c0.x)/nc.x;
        a1 = ns.x/nc.x;

        a2 = s0.y - c0.y - a0*nc.y;
        a3 = a1*nc.y - ns.y;

        if ( fabs(a3) > 1.0e-5 ) //parallel check
        {
          ts = a2/a3;
          tc = a0 + a1*ts;

          if ( (ts < (1.0-cut_tol)) && (ts > (0.0+cut_tol)) )
          {
            //printf("Cutting segment (%f,%f)->(%f,%f)\n", s0.x,s0.y,s1.x,s1.y);
            //========================= Erase intersected segment
            surface_mesh->lines.erase(surface_mesh->lines.begin()+ell);

            //========================= Build new lines
            Vertex new_vert = s0 + ns*ts;
            surface_mesh->vertices.push_back(new_vert);
            int new_vert_index = surface_mesh->vertices.size()-1;
            inserted_verts.push_back(new_vert_index);

            //printf("Inserted vertex %d %f,%f\n",new_vert_index,new_vert.x,new_vert.y);

            Line new_line_A, new_line_B;

            new_line_A.vertices[0] = s0;
            new_line_A.vertices[1] = new_vert;

            new_line_A.v_index[0] = s0_index;
            new_line_A.v_index[1] = new_vert_index;

            new_line_B.vertices[0] = new_vert;
            new_line_B.vertices[1] = s1;

            new_line_B.v_index[0] = new_vert_index;
            new_line_B.v_index[1] = s1_index;

            surface_mesh->lines.push_back(new_line_A);
            surface_mesh->lines.push_back(new_line_B);

            found_intersection = true;
            break;

          } //if in tol
        } //if not parallel
      } //if horizontal line

    }//for segment

  }//while

  //================================================== Sort inserted vertices
  Vertex LL(-1.0e6,-1.0e6);
  for (int i=0; i<(inserted_verts.size()-1); i++)
  {


    for (int j=0; j<(inserted_verts.size()-1-i); j++)
    {
      Vertex vj = surface_mesh->vertices[inserted_verts[j]] - LL;
      double vjd = vj.Norm();

      Vertex vjp1 = surface_mesh->vertices[inserted_verts[j+1]] - LL;
      double vjp1d = vjp1.Norm();

      if (vjd > vjp1d)
      {
//        printf("Swapping %d and %d\n", inserted_verts[j], inserted_verts[j+1]);
//        printf("yj=%f, yjp1=%f\n",vj.y,vjp1.y);
        int temp = inserted_verts[j+1];
        inserted_verts[j+1] = inserted_verts[j];
        inserted_verts[j]=temp;

      }
    }
  }

//  for (int i=0; i<(inserted_verts.size()); i++)
//  {
//    printf("Sorted vertex %d, index=%d x=%f\n",i,inserted_verts[i],surface_mesh->vertices[inserted_verts[i]].y);
//  }

  //================================================== Fill in actual cut line
  for (unsigned v=0; v<(inserted_verts.size()-1); v++)
  {
    Line new_line;
    int vi_index = inserted_verts[v];
    int vf_index = inserted_verts[v+1];

    new_line.vertices[0] = surface_mesh->vertices[vi_index];
    new_line.vertices[1] = surface_mesh->vertices[vf_index];

    new_line.v_index[0] = vi_index;
    new_line.v_index[1] = vf_index;

    surface_mesh->lines.push_back(new_line);
  }

  //================================================== Clean up
  inserted_verts.clear();
}
