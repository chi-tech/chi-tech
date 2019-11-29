#ifndef _chi_meshhandler_h
#define _chi_meshhandler_h

#include<stdio.h>
#include <vector>

#include"../chi_mesh.h"
#include "../LineMesh/chi_linemesh.h"


/**Object for containing all mesh related variables.*/
class chi_mesh::MeshHandler
{
public:
  typedef std::vector<chi_mesh::SurfaceMesh*> SurfaceMeshCollection;
  std::vector<SurfaceMeshCollection*>      surface_mesh_collections;
  std::vector<chi_mesh::SurfaceMesh*>      surface_mesh_stack;
  std::vector<chi_mesh::Region*>           region_stack;
  std::vector<chi_mesh::LineMesh*>         linemesh_stack;
  std::vector<chi_mesh::LogicalVolume*>    logicvolume_stack;
  std::vector<chi_mesh::FieldFunctionInterpolation*> ffinterpolation_stack;

  std::vector<chi_mesh::EdgeLoopCollection*> edge_loop_collections;

  chi_mesh::SurfaceMesher* surface_mesher;
  chi_mesh::VolumeMesher*  volume_mesher;

  std::vector<chi_mesh::Cell*> pCells;

  std::vector<chi_mesh::CELL_SET*> cell_sets;

public:
  chi_mesh::MeshContinuum* GetGrid();
};

#endif