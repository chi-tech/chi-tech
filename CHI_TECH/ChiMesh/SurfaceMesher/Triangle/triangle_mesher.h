#ifndef _triangle_mesher_h
#define _triangle_mesher_h

#include "../surfacemesher.h"

//######################################################### Class def
/**Triangulation of a mesh.
 *


 */
class chi_mesh::SurfaceMesherTriangle : public chi_mesh::SurfaceMesher
{
public:
  bool get_auto_min;
  double area_constraint;



  double cut_tol;
public:
  //00
  SurfaceMesherTriangle() : SurfaceMesher(SurfaceMesherType::Delaunay)
  {
    get_auto_min = true;
    area_constraint = 1.0;
    partitioning_x = 1;
    partitioning_y = 1;

    cut_tol = 0.005;

    export_loadbalance = false;
  }
  //01
  void Execute();
  //02 Utils
  void CalculateMinimumArea(chi_mesh::Boundary* boundary);
  void ReorderCells(chi_mesh::SurfaceMesh* surface_mesh,
                    chi_mesh::SurfaceMesh* out_surface_mesh);
  //03
  void MeshBoundary(chi_mesh::Boundary* boundary);
  //04
  void AddCutLines(chi_mesh::Boundary* boundary);
  //04a
  void IntersectSegments(chi_mesh::Line cut_line,
                         chi_mesh::SurfaceMesh* surface_mesh);
};


#endif