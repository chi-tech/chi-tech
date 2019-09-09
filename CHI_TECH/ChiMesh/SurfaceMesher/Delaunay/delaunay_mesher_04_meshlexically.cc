#include "delaunay_mesher.h"
#include <unistd.h>

//##############################################################
/**Performs a lexicographical meshing of a patch.*/
void chi_mesh::SurfaceMesherDelaunay::MeshLexically(
  chi_mesh::SurfaceMesherDelaunay::Patch *raw_patch,
  chi_mesh::SurfaceMesherDelaunay::DelaunayPatch &del_patch)
{
  //================================================== Connect parent patch
  del_patch.raw_patch = raw_patch;

  //================================================== Generate simplices
  del_patch.GenerateSimplices();

  //================================================== Assign normal
  chi_mesh::Face first_face = raw_patch->faces.front();
  del_patch.patch_normal = first_face.geometric_normal;

  //================================================== Project to 2D
  del_patch.ProjectTo2D();

  //================================================== Lexicographicall sort
  //                                                   vertices
  del_patch.LexicographicallySortVertices();

  //================================================== Create first triangle
  del_patch.CreateFirstTriangle();

  //================================================== Attach unused vertices
  while (del_patch.unused_vs.size()>0)
  {
    int v_index = del_patch.unused_vs.back();

    del_patch.AttachUnusedVertex(v_index);

    while (del_patch.ConvexifyHull())
    {

    }

//    std::vector<int>::iterator unused;
//    for (unused = del_patch.unused_vs.begin();
//         unused != del_patch.unused_vs.end();
//         unused++)
//    {
//      int index = *unused;
//      Vertex v = del_patch.Pstar.at(index);
//      printf("Unused vertex %d %+.4f %+.4f %+.4f\n",index,v.x,v.y,v.z);
//    }

    std::vector<Edge>::iterator curedge;
    for (curedge = del_patch.open_edges.begin();
         curedge!=del_patch.open_edges.end();
         curedge++)
    {
      printf("Open edge %d->%d    ", curedge->v_index[0],curedge->v_index[1]);
      printf("(%+.3f,%+.3f)->(%+.3f,%+.3f)\n",curedge->vertices[0].x,
                                              curedge->vertices[0].y,
                                              curedge->vertices[1].x,
                                              curedge->vertices[1].y);
    }
    printf("Number of unused vertices = %d\n",del_patch.unused_vs.size());

  }


  del_patch.ExportAsObj("CHI_TEST/Zmesh.obj");

}