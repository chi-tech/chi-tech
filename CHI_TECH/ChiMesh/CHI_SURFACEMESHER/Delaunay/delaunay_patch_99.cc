#include "delaunay_mesher.h"

//###################################################################
/**Exports a delaunay patch to file. This is meant to be used for
 * debug purposes only.*/
void chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::ExportAsObj(
  const char *fileName)
{
  if (this->triangles.size()==0)
  {
    std::cout << "Cannot export empty DelaunayPatch\n";
    return;
  }
  FILE* outputFile = fopen(fileName,"w");
  if (outputFile==NULL)
  {
    printf("Error creating file %s!\n",fileName);
    return;
  }

  fprintf(outputFile,"# Exported mesh file from tringulation script\n");
  fprintf(outputFile,"o %s\n","Chitech Delaunay Patch");

  std::vector<chi_mesh::Vertex>::iterator cur_v;
  for (cur_v = this->Pstar.begin();
       cur_v != this->Pstar.end();
       cur_v++)
  {
    fprintf(outputFile,"v %9.6f %9.6f %9.6f\n",cur_v->x,cur_v->y,cur_v->z);
  }

  chi_mesh::Face first_face = this->triangles.front();
  fprintf(outputFile,"vn %.4f %.4f %.4f\n", first_face.geometric_normal.x,
          first_face.geometric_normal.y,
          first_face.geometric_normal.z);
  fprintf(outputFile,"s off\n");

  std::vector<chi_mesh::Face>::iterator cur_face;
  for (cur_face = this->triangles.begin();
       cur_face != this->triangles.end();
       cur_face++)
  {
    fprintf(outputFile,"f %d//1 %d//1 %d//1\n",cur_face->v_index[0]+1,
            cur_face->v_index[1]+1,
            cur_face->v_index[2]+1);
  }

  fclose(outputFile);
  printf("Exported mesh to %s\n",fileName);
}