#include"chi_mesh.h"
#include"SurfaceMesh/chi_surfacemesh.h"
#include"Boundary/chi_boundary.h"


chi_mesh::Boundary* chi_mesh::AssignSurfaceToBoundary(chi_mesh::SurfaceMesh *surface)
{
  chi_mesh::Boundary* new_boundary = new chi_mesh::Boundary;

  /*//======================================================= Copy vertices
  std::vector<chi_mesh::Vertex>::iterator curVertex;
  for(curVertex = surface->vertices.begin(); curVertex != surface->vertices.end(); curVertex++)
  {
    chi_mesh::Vertex new_vertex = *curVertex;
    new_boundary->vertices.push_back(new_vertex);
  }

  //======================================================= Copy texture vertices
  for(curVertex = surface->tex_vertices.begin(); curVertex != surface->tex_vertices.end(); curVertex++)
  {
    chi_mesh::Vertex new_vertex = *curVertex;
    new_boundary->tex_vertices.push_back(new_vertex);
  }

  //======================================================= Copy normals
  std::vector<chi_mesh::Normal>::iterator curNormal;
  for(curVertex = surface->normals.begin(); curVertex != surface->normals.end(); curVertex++)
  {
    chi_mesh::Normal new_normal = *curNormal;
    new_boundary->normals.push_back(new_normal);
  }

  //======================================================= Copy faces
  std::vector<chi_mesh::Face>::iterator curFace;
  for (curFace = surface->faces.begin(); curFace != surface->faces.end(); curFace++)
  {
    chi_mesh::Face new_face = *curFace;
    new_boundary->faces.push_back(new_face);
  }*/

  return new_boundary;
}