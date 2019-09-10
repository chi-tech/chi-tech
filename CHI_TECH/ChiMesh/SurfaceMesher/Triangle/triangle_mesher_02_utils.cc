#include "triangle_mesher.h"
#include "../../MeshContinuum/chi_meshcontinuum.h"
#include "../../Boundary/chi_boundary.h"

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**Executes the triangle mesher*/
void chi_mesh::SurfaceMesherTriangle::CalculateMinimumArea(chi_mesh::Boundary *boundary)
{
  chi_mesh::SurfaceMesh* surf_mesh = boundary->initial_mesh_continuum.surface_mesh;

  //================================================== Compute centroid
  chi_mesh::Vertex centroid;
  std::vector<chi_mesh::Vertex>::iterator curvert;
  for (curvert = surf_mesh->vertices.begin();
       curvert != surf_mesh->vertices.end();
       curvert++)
  {
    centroid = centroid + (*curvert);
  }
  centroid = centroid/surf_mesh->vertices.size();

  //================================================== Calculate max distance
  //                                                   from centroid
  double maxD = 0.0;
  for (curvert = surf_mesh->vertices.begin();
       curvert != surf_mesh->vertices.end();
       curvert++)
  {
    double d = ((*curvert) - centroid).Norm();
    if (d>maxD) {maxD = d;}
  }

  //================================================== Compute area constraint
  this->area_constraint = (maxD/20.0)*(maxD/20.0);
}


//###################################################################
/**Reorders the faces of a created surface mesh.*/
void chi_mesh::SurfaceMesherTriangle::
  ReorderCells(chi_mesh::SurfaceMesh *surface_mesh,
               chi_mesh::SurfaceMesh* out_surface_mesh)
{
  chi_log.Log(LOG_ALLVERBOSE_2) << "Reordering triangle order";
  //============================================= Create the graph
  CHI_UD_GRAPH ugraph;
  std::vector<int> mapping;

  //============================================= Add a vertex for each face
  //                                              and init mapping
  for (int f=0; f<surface_mesh->faces.size(); f++)
  {
    boost::add_vertex(ugraph);
    mapping.push_back(f);
  }

  //============================================= Connect the vertices
  for (int f=0; f<surface_mesh->faces.size(); f++)
  {
    chi_mesh::Face* curface = &surface_mesh->faces[f];
    for (int e=0; e<3; e++)
    {
      if (curface->e_index[e][2]>=0)
      {
        boost::add_edge(f,curface->e_index[e][2],ugraph);
      }
    }
  }

  //============================================= Perform reordering
  chi_graph::CuthillMckee(ugraph,&mapping);

  //============================================= InitializeAlphaElements the out mesh
  //=================================== Copy vertices
  for (int v=0; v<surface_mesh->vertices.size(); v++)
  {
    out_surface_mesh->vertices.push_back(surface_mesh->vertices[v]);
  }
  //=================================== Copy faces
  for (int f=0; f<mapping.size(); f++)
  {
    int out_index = mapping[f];
    out_surface_mesh->faces.push_back(surface_mesh->faces[out_index]);
  }
//  out_surface_mesh = surface_mesh;
}


