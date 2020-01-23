#ifndef _chi_meshcontinuum_h
#define _chi_meshcontinuum_h

#include "../chi_mesh.h"
#include <boost/graph/adjacency_list.hpp>
#include "../../ChiGraph/chi_graph.h"
#include "../Cell/cell.h"



//######################################################### Class Definition
/**Stores the relevant information for completely defining a computational
 * domain. */
class chi_mesh::MeshContinuum
{
public:
  //################################################## LocalCells Class Definition
  /**Stores references to global cells to enable an iterator.*/
  class LocalCells
  {
  public:
    std::vector<int>& local_cell_ind;
    std::vector<chi_mesh::Cell*>& cell_references;

    LocalCells(std::vector<int>& in_local_cell_ind,
               std::vector<chi_mesh::Cell*>& in_cell_references) :
      local_cell_ind(in_local_cell_ind),
      cell_references(in_cell_references) {}

    //##################################### iterator Class Definition
    /**Internal iterator class.*/
    class iterator
    {
    public:
      LocalCells& ref_block;
      size_t      ref_element;

      iterator(LocalCells& in_block, size_t i) :
        ref_block(in_block),
        ref_element(i) {}

      iterator operator++(        ) {iterator i = *this; ref_element++; return i;}
      iterator operator++(int junk) {ref_element++; return *this;}
      chi_mesh::Cell& operator*()
      {
        return *(ref_block.cell_references[
                   ref_block.local_cell_ind[ref_element]]);
      }
      chi_mesh::Cell* operator->()
      {
        return ref_block.cell_references[
                 ref_block.local_cell_ind[ref_element]];}
      bool operator==(const iterator& rhs)
      {
        return ref_element == rhs.ref_element;
      }
      bool operator!=(const iterator& rhs)
      {
        return ref_element != rhs.ref_element;
      }
    };

    iterator begin() {return iterator(*this,0);}

    iterator end(){return iterator(*this,local_cell_ind.size());}
  };
  //--------------------------------------------------

  std::vector<chi_mesh::Node*>   nodes;
  std::vector<chi_mesh::Cell*>   cells;
  LocalCells                     local_cells;
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
  MeshContinuum() : local_cells(local_cell_glob_indices,cells)
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

  int  FindAssociatedFace(chi_mesh::CellFace& cur_face,int adj_cell_g_index,bool verbose=false);
  void FindAssociatedVertices(chi_mesh::CellFace& cur_face,
                              int adj_cell_g_index, int associated_face,
                              std::vector<int>& dof_mapping);
  void PopulateUniqueBoundaries(std::set<int>& bndries);


};




#endif