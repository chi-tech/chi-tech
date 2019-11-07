#ifndef _chi_meshcontinuum_h
#define _chi_meshcontinuum_h

#include "../chi_mesh.h"
#include <boost/graph/adjacency_list.hpp>
#include "../../ChiGraph/chi_graph.h"
#include "../Cell/cell.h"


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

private:
  bool                           face_histogram_available;

  //Pair.first is the max dofs-per-face for the category and Pair.second
  //is the number of faces in this category
  std::vector<std::pair<size_t,size_t>> face_categories;

public:
  MeshContinuum()
  {
    this->surface_mesh = nullptr;
    this->line_mesh    = nullptr;
    face_histogram_available = false;
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
  void BuildFaceHistogramInfo(double master_tolerance=1.2, double slave_tolerance=1.1);
  size_t NumberOfFaceHistogramBins();
  size_t MapFaceHistogramBins(size_t num_face_dofs);
  size_t GetFaceHistogramBinDOFSize(size_t category);
  bool IsCellLocal(int cell_global_index=-1);
  bool IsCellBndry(int cell_global_index = 0);

//  int  FindAssociatedFace(chi_mesh::PolyFace* cur_face,int adj_cell_g_index,bool verbose_info=false);
  int  FindAssociatedFace(chi_mesh::CellFace& cur_face,int adj_cell_g_index,bool verbose=false);
//  int  FindAssociatedEdge(int* edgeinfo,int adj_cell_g_index,bool verbose_info=false);
//  void FindAssociatedVertices(chi_mesh::PolyFace* cur_face,
//                              int adj_cell_g_index, int associated_face,
//                              std::vector<int>& dof_mapping);
  void FindAssociatedVertices(chi_mesh::CellFace& cur_face,
                              int adj_cell_g_index, int associated_face,
                              std::vector<int>& dof_mapping);


};


#endif