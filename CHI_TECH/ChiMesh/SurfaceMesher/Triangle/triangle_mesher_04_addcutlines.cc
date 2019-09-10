#include "triangle_mesher.h"
#include "../../MeshContinuum/chi_meshcontinuum.h"
#include "../../Boundary/chi_boundary.h"

//###################################################################
/**Adds cut lines to the surface mesh of a boundary.*/
void chi_mesh::SurfaceMesherTriangle::AddCutLines(chi_mesh::Boundary *boundary)
{
  //================================================== Get reference to surface
  chi_mesh::SurfaceMesh* surface_mesh =
    boundary->initial_mesh_continuum.surface_mesh;

  //================================================== Check correct xcut amount
  if (partitioning_x>1)
  {
    if (xcuts.size()!= (partitioning_x-1))
    {
      std::cerr << "Incorrect amount of x-cuts (";
      std::cerr << xcuts.size();
      std::cerr << ") for the given x partitioning intervals (";
      std::cerr << partitioning_x << ")\n";
      exit(EXIT_FAILURE);
    }
  }

  //================================================== Check correct ycut amount
  if (partitioning_y>1)
  {
    if (ycuts.size()!= (partitioning_y-1))
    {
      std::cerr << "Incorrect amount of y-cuts (";
      std::cerr << ycuts.size();
      std::cerr << ") for the given y partitioning intervals (";
      std::cerr << partitioning_y << ")\n";
      exit(EXIT_FAILURE);
    }
  }

  //================================================== Find minimums & maximums
  double max_x = -1.0e10;
  double max_y = -1.0e10;

  double min_x = 1.0e10;
  double min_y = 1.0e10;
  for (unsigned v=0; v<surface_mesh->vertices.size(); v++)
  {
    if (surface_mesh->vertices[v].x > max_x)
    {max_x = surface_mesh->vertices[v].x;}

    if (surface_mesh->vertices[v].y > max_y)
    {max_y = surface_mesh->vertices[v].y;}

    if (surface_mesh->vertices[v].x < min_x)
    {min_x = surface_mesh->vertices[v].x;}

    if (surface_mesh->vertices[v].y < min_y)
    {min_y = surface_mesh->vertices[v].y;}
  }

  //================================================== Check legal cuts-x
  for (unsigned i=0; i < xcuts.size(); i++)
  {
    if ( (xcuts[i] < min_x) || (xcuts[i] > max_x) )
    {
      fprintf(stderr, "X-cut %d, at x=%f, is outside mesh "
                      "(x in [%f,%f])\n",
              i, xcuts[i],min_x,max_x);
    }
  }

  //================================================== Check legal cuts-y
  for (unsigned i=0; i < ycuts.size(); i++)
  {
    if ( (ycuts[i] < min_y) || (ycuts[i] > max_y) )
    {
      fprintf(stderr, "Y-cut %d, at y=%f, is outside mesh "
                      "(y in [%f,%f])\n",
              i, ycuts[i],min_y,max_y);
    }
  }

  //================================================== Compute xcuts segment
  //                                                   intersections
  for (unsigned i=0; i<xcuts.size(); i++)
  {
    Line cut_line;
    Vertex vi(xcuts[i],min_y);
    Vertex vf(xcuts[i],max_y);

    cut_line.vertices[0] = vi;
    cut_line.vertices[1] = vf;

    //printf("Cutline (%f,%f)->(%f,%f)\n", vi.x,vi.y,vf.x,vf.y);

    IntersectSegments(cut_line, surface_mesh);
  }


  for (unsigned i=0; i<ycuts.size(); i++)
  {
    Line cut_line;
    Vertex vi(min_x,ycuts[i]);
    Vertex vf(max_x,ycuts[i]);

    cut_line.vertices[0] = vi;
    cut_line.vertices[1] = vf;

    //printf("Cutline (%f,%f)->(%f,%f)\n", vi.x,vi.y,vf.x,vf.y);

    IntersectSegments(cut_line, surface_mesh);
  }

}