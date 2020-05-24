#ifndef _chi_meshcontinuum_h
#define _chi_meshcontinuum_h

#include "../chi_mesh.h"
#include <boost/graph/adjacency_list.hpp>
#include "../../ChiGraph/chi_graph.h"
#include "ChiMesh/Cell/cell.h"

#include <chi_mpi.h>

//######################################################### Class Definition
/**Stores the relevant information for completely defining a computational
 * domain. */
class chi_mesh::MeshContinuum
{
public:
  //##################################################
  /**Stores references to global cells to enable an iterator.*/
  class LocalCells
  {
  public:
    std::vector<int>& local_cell_ind;

    std::vector<chi_mesh::Cell*> native_cells;
    std::vector<chi_mesh::Cell*> foreign_cells;

    /**Constructor.*/
    LocalCells(std::vector<int>& in_local_cell_ind) :
      local_cell_ind(in_local_cell_ind)
      {}

    ~LocalCells()
    {
      for (auto cell : native_cells) delete cell;
      for (auto cell : foreign_cells) delete cell;
    }

    chi_mesh::Cell& operator[](int cell_local_index);


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
      { return *(ref_block.native_cells[ref_element]); }
      chi_mesh::Cell* operator->()
      { return ref_block.native_cells[ref_element]; }
      bool operator==(const iterator& rhs)
      { return ref_element == rhs.ref_element; }
      bool operator!=(const iterator& rhs)
      { return ref_element != rhs.ref_element; }
    };

    iterator begin() {return {*this,0};}

    iterator end(){return {*this,local_cell_ind.size()};}

    size_t size() {return local_cell_ind.size();}
  };
  //--------------------------------------------------

  //##################################################
  /**Handles all global index queries.*/
  class GlobalCellHandler
  {
  private:
    LocalCells& local_cells;
    std::set<std::pair<int,int>> global_cell_native_index_set;
    std::set<std::pair<int,int>> global_cell_foreign_index_set;


  public:
    explicit GlobalCellHandler(LocalCells& local_cell_object_ref) :
     local_cells(local_cell_object_ref)
    {}

    void push_back(chi_mesh::Cell* new_cell);
    chi_mesh::Cell* &operator[](int cell_global_index);


  };

  std::vector<chi_mesh::Node*>   vertices;
  LocalCells                     local_cells;
  GlobalCellHandler              cells;
  chi_mesh::SurfaceMesh*         surface_mesh;
  chi_mesh::LineMesh*            line_mesh;
  std::vector<int>               local_cell_glob_indices;
  std::vector<int>               boundary_cell_indices;




private:
  bool                           face_histogram_available = false;
  bool                           communicators_available  = false;

  //Pair.first is the max dofs-per-face for the category and Pair.second
  //is the number of faces in this category
  std::vector<std::pair<size_t,size_t>> face_categories;

  ChiMPICommunicatorSet commicator_set;

public:
  MeshContinuum() :
    local_cells(local_cell_glob_indices),
    cells(local_cells)
  {
    surface_mesh = nullptr;
    line_mesh    = nullptr;
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

  void FindAssociatedVertices(chi_mesh::CellFace& cur_face,
                              std::vector<int>& dof_mapping);

  chi_mesh::Vector3 ComputeCentroidFromListOfNodes(const std::vector<int>& list);

  void CommunicatePartitionNeighborCells(
    std::vector<chi_mesh::Cell*>& neighbor_cells);

  ChiMPICommunicatorSet& GetCommunicator();

  size_t GetGlobalNumberOfCells();
};




#endif