#include "delaunay_mesher.h"
#include <unistd.h>

//##############################################################
/**Performs a lexicographical meshing of a patch.*/
void chi_mesh::SurfaceMesherDelaunay::MeshDelaunay(
  chi_mesh::SurfaceMesherDelaunay::Patch *raw_patch,
  chi_mesh::SurfaceMesherDelaunay::DelaunayPatch &del_patch)
{
  //================================================== Connect parent patch
  del_patch.raw_patch = raw_patch;

  //================================================== Copy vertices
  for (unsigned v=0;v<raw_patch->vertices.size(); v++)
  {
    del_patch.vertices.push_back(raw_patch->vertices[v]);
  }

  //================================================== Copy faces
  for (unsigned f=0;f<raw_patch->faces.size(); f++)
  {
    del_patch.triangles.push_back(raw_patch->faces[f]);
  }

  //================================================== Generate simplices
  del_patch.GenerateSimplices();

  //================================================== Assign normal
  chi_mesh::Face first_face = raw_patch->faces.front();
  del_patch.patch_normal = first_face.geometric_normal;

  //================================================== Project to 2D
  del_patch.ProjectTo2D();

  //================================================== Restore Delaunay condition
  del_patch.RestoreDelaunay();

  //================================================== Refine the mesh
  int encroach_iter=0;
  while (del_patch.RemoveEncroachedEdges())
  {
    del_patch.RestoreDelaunay();
    encroach_iter++;
    //if (encroach_iter>4) break;
  }



  bool dont_stop = true;
  for (int k=0; k<300; k++)
  {

    printf("Refinement iteration %d\n",k);
    dont_stop = del_patch.RemoveBadQualityTriangles();
    if (dont_stop)
    {
      del_patch.RestoreDelaunay();
      //del_patch.DumpToScilab("ZScilabMesh.sce");
      while (del_patch.RemoveEncroachedEdges())
      {

        del_patch.RestoreDelaunay();
      }
    }

    printf("Refine Edges\n");


    dont_stop = del_patch.RefineEdges();
    if (dont_stop)
    {
      del_patch.RestoreDelaunay();
      //del_patch.DumpToScilab("ZScilabMesh.sce");
      while (del_patch.RemoveEncroachedEdges())
      {

        del_patch.RestoreDelaunay();

      }
    }

    if (!dont_stop) break;
  }

  //del_patch.ExportAsObj("CHI_TEST/Zmesh.obj");
  del_patch.DumpToScilab("ZScilabMesh.sce");

}