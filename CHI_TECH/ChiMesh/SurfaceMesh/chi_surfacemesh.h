#ifndef _chi_surfacemesh_h
#define _chi_surfacemesh_h

#include<stdio.h>
#include <vector>

#include"../chi_mesh.h"



/**Generic surface mesh class.*/
class chi_mesh::SurfaceMesh
{
public:
  std::vector<chi_mesh::Vertex>   vertices;
  std::vector<chi_mesh::Vertex>   tex_vertices;
  std::vector<chi_mesh::Normal>   normals;
  std::vector<chi_mesh::Face>     faces;
  std::vector<chi_mesh::Edge>     lines;
  std::vector<chi_mesh::PolyFace*> poly_faces;

public:
  //00
        SurfaceMesh();
  friend std::ostream& operator<<(std::ostream& os,  SurfaceMesh& dt);
  //01
  int   ImportFromOBJFile(const char* fileName,bool as_poly);
  int   ImportFromTriangleFiles(const char* fileName, bool as_poly);
  void  UpdateInternalConnectivity();
  void  ExportToOBJFile(const char* fileName);
  void  ExportToPolyFile(const char* fileName);

  //02
  bool  CheckNegativeSense(double x, double y, double z);

  //03
  EdgeLoopCollection*  GetEdgeLoops();
  EdgeLoopCollection*  GetCoLinearEdges(EdgeLoopCollection* in_loop);
  //03b
  EdgeLoopCollection*  GetEdgeLoopsPoly();

  //04
  void  SplitByPatch(std::vector<chi_mesh::SurfaceMesh*>& patches);

  //05
  void  ExtractOpenEdgesToObj(const char* fileName);

  //06
  void  CheckCyclicDependencies(int num_angles);
  void  GetMeshStats();

};

#endif