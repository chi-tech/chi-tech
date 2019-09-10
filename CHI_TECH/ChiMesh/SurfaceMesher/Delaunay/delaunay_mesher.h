#ifndef _delaunay_mesher_h
#define _delaunay_mesher_h

#include <iostream>
#include<vector>
#include "../surfacemesher.h"
#include "../../SurfaceMesh/chi_surfacemesh.h"




//######################################################### Class def
/**Delaunay triangulation of a mesh.
 *
## Step 1 Create contexts
Create a Delaunay context for each boundary. A delaunay context
 is a data structure that will hold the delaunay mesh that is to be
 generated. Line meshes are ignored but still parsed.

 \code
   struct DelaunayMeshContext
  {
    chi_mesh::Boundary*    context_boundary;
    chi_mesh::SurfaceMesh* context_surface_mesh;
    chi_mesh::LineMesh*    context_line_mesh;

    std::vector<chi_mesh::SurfaceMesh*> patches;
    std::vector<chi_mesh::SurfaceMesh*> delaunay_patches;

    std::vector<chi_mesh::Interface*> interfaces;
  };
 \endcode

 ## Step 2 Mesh each context
 For each surface context an attempt is made to split it by patch. This is
 necessary because the 2D algorithm cannot handle surfaces that are not
 co-planar.

 */
class chi_mesh::SurfaceMesherDelaunay : public chi_mesh::SurfaceMesher
{
public:
  /**This context assists in three situations. The first is where there are
 * multiple boundaries between regions that requires conformal meshing (connected
 * by edges). The second
 * is where a single boundary comprises a number of non-planar non-contigous
 * surfaces where each non planar surface need to be conformal as well
 * (through edge connections). The third type is interfaces between regions
 * which can be adjoined through edges and/or surfaces. */
  struct DelaunayMeshContext
  {
    chi_mesh::Boundary*    context_boundary;
    chi_mesh::SurfaceMesh* context_surface_mesh;
    chi_mesh::LineMesh*    context_line_mesh;

    std::vector<chi_mesh::SurfaceMesh*> patches;
    std::vector<chi_mesh::SurfaceMesh*> delaunay_patches;

    std::vector<chi_mesh::Interface*> interfaces;

    chi_mesh::SurfaceMesh* remeshed_surface_mesh;
  };

  typedef std::vector<DelaunayMeshContext*> RegionDelaunyContexts;
  typedef chi_mesh::SurfaceMesh Patch;
  typedef chi_mesh::Face Tri;
  //typedef chi_mesh::SurfaceMesh DelaunayPatch;
  class DelaunayPatch;

  double      tolerance;
  std::vector<RegionDelaunyContexts*> region_contexts;
public:
  //00
  SurfaceMesherDelaunay();
  //02
  void Execute();
  //03
  void MeshRegion(RegionDelaunyContexts& region_context);
  //04
  void MeshLexically(Patch* raw_patch, DelaunayPatch& del_patch);
  //05
  void MeshDelaunay(Patch* raw_patch, DelaunayPatch& del_patch);
};

//#########################################################
/**Handles data and methods specific to Delaunay triangulation.*/
class chi_mesh::SurfaceMesherDelaunay::DelaunayPatch :
  public chi_mesh::SurfaceMesh
{
public:
  Patch*              raw_patch;
  std::vector<int*>   lexi_map;
  EdgeLoopCollection  simplices;
  chi_mesh::Normal    patch_normal;
  chi_mesh::Vector    hat_i,hat_j,hat_k;
  chi_mesh::Vector    centroid;
  std::vector<chi_mesh::Vertex> Pstar;
  std::vector<int>    unused_vs;

  std::vector<Tri>    triangles;
  std::vector<Edge>   open_edges;
  //std::vector<Edge>   non_loc_del_edges;

  double tolerance;
  double absoluteMinumumSize;
  double quality_factor1;
  bool get_auto_min;
public:
  //00 utils
  DelaunayPatch(double tol)
  {
    tolerance = tol;
    absoluteMinumumSize = 0.1;
    quality_factor1 = 1.0;
    get_auto_min = true;
  }
  double Orient2D(Vertex a, Vertex b, Vertex c);
  bool   CheckCrossSimplices(Vertex v0, Vertex v1);
  double InCircle(Vertex a,Vertex b,Vertex c,Vertex d);
  void   FindNonLocDelEdge(std::vector<Edge>& non_loc_del_edges);
  void   DumpToScilab(const char* file_name, bool verbose=false);
  bool   CheckEdgeCanBeSplit(unsigned i_vindex,unsigned f_vindex);
  bool   CheckEdgeEncroached(unsigned i_vindex,
                             unsigned f_vindex,
                             unsigned other_vindex);
  bool   CheckEdgeEncroached(unsigned i_vindex,
                             unsigned f_vindex,
                             Vertex v2);
  void FindCircumCircle(Vertex a, Vertex b, Vertex c,
                        Vertex& center, double& radius);
  //01
  void GenerateSimplices(bool add_vertices=false);
  //02
  void ProjectTo2D();
  //03
  void LexicographicallySortVertices();
  //04
  void CreateFirstTriangle();
  //05
  void AttachUnusedVertex(int v_index);
  //06
  bool ConvexifyHull();
  //07
  void RestoreDelaunay();
  void EdgeFlip(Edge& edge_to_flip);
  bool RemoveEncroachedEdges();
  bool InsertVertex(Vertex vc_pstar, unsigned seed_search_triangle_index);
  bool RemoveBadQualityTriangles();
  bool RefineEdges();

  //99
  void ExportAsObj(const char *fileName);
};

#endif