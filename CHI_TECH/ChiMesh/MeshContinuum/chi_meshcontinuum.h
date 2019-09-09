#ifndef _chi_meshcontinuum_h
#define _chi_meshcontinuum_h

#include "../chi_mesh.h"
#include <boost/graph/adjacency_list.hpp>
#include "../../ChiGraph/chi_graph.h"


//######################################################### Class Definition
class chi_mesh::MeshContinuum
{
public:
  std::vector<chi_mesh::Node*>   nodes;
  std::vector<chi_mesh::Cell*>   cells;
  chi_mesh::SurfaceMesh*         surface_mesh;
  chi_mesh::LineMesh*            line_mesh;
  std::vector<int>               local_cell_glob_indices;
  std::vector<int>               glob_cell_local_indices;
  std::vector<int>               boundary_cell_indices;

  //CHI_UD_GRAPH                   grid_graph;

  MeshContinuum()
  {
    this->surface_mesh = nullptr;
    this->line_mesh    = nullptr;

  }

  //01
  void ExportCellsToPython(const char* fileName,
                           bool surface_only=true,
                           std::vector<int>* cell_flags = nullptr,
                           int options = 0);
  void ExportCellsToObj(const char* fileName,
                           bool per_material=false,
                           int options = 0);
  void ExportCellsToVTK(const char* baseName);

  //02
  void ConnectGrid();
  bool IsCellLocal(int cell_global_index=-1);
  bool IsCellBndry(int cell_global_index = 0);

  int  FindAssociatedFace(chi_mesh::PolyFace* cur_face,int adj_cell_g_index,bool verbose=false);
  int  FindAssociatedEdge(int* edgeinfo,int adj_cell_g_index,bool verbose=false);
  void FindAssociatedVertices(chi_mesh::PolyFace* cur_face,
                              int adj_cell_g_index, int associated_face,
                              std::vector<int>& dof_mapping);


};


#endif