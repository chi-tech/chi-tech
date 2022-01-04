#ifndef CHI_MESHCONTINUUM_H_
#define CHI_MESHCONTINUUM_H_

#include "../chi_mesh.h"
#include "chi_meshcontinuum_localcellhandler.h"
#include "chi_meshcontinuum_globalcellhandler.h"
#include "chi_meshcontinuum_vertexhandler.h"

#include "chi_mpi.h"

//######################################################### Class Definition
/**Stores the relevant information for completely defining a computational
 * domain. */
class chi_mesh::MeshContinuum
{
private:
  std::vector<chi_mesh::Cell*> native_cells;  ///< Actual native cells
  std::vector<chi_mesh::Cell*> foreign_cells; ///< Locally stored ghosts

  std::map<uint64_t,uint64_t> global_cell_id_to_native_id_map;
  std::map<uint64_t,uint64_t> global_cell_id_to_foreign_id_map;

  uint64_t global_vertex_count=0;

public:
  VertexHandler                  vertices;
  LocalCellHandler               local_cells;
  GlobalCellHandler              cells;
  chi_mesh::SurfaceMesh*         surface_mesh;
  chi_mesh::LineMesh*            line_mesh;
  std::vector<uint64_t>          local_cell_glob_indices;

private:
  bool                           face_histogram_available = false;
  bool                           communicators_available  = false;

  //Pair.first is the max dofs-per-face for the category and Pair.second
  //is the number of faces in this category
  std::vector<std::pair<size_t,size_t>> face_categories;

  ChiMPICommunicatorSet commicator_set;

public:
  MeshContinuum() :
    local_cells(native_cells, foreign_cells),
    cells(local_cell_glob_indices,
          native_cells,
          foreign_cells,
          global_cell_id_to_native_id_map,
          global_cell_id_to_foreign_id_map)
  {
    surface_mesh = nullptr;
    line_mesh    = nullptr;
  }

  void SetGlobalVertexCount(const uint64_t count) {global_vertex_count = count;}
  uint64_t GetGlobalVertexCount() const {return global_vertex_count;}

  static
  std::shared_ptr<MeshContinuum> New()
  { return std::shared_ptr<MeshContinuum>(new MeshContinuum());}

  /**Method to be called if cells and nodes have been transferred
   * to another grid.*/
  void ClearCellReferences()
  {
    native_cells.clear();
    foreign_cells.clear();
    global_cell_id_to_native_id_map.clear();
    global_cell_id_to_foreign_id_map.clear();
  }

  //01
  void ExportCellsToPython(const char* fileName,
                           bool surface_only=true,
                           std::vector<int>* cell_flags = nullptr,
                           int options = 0) const;
  void ExportCellsToObj(const char* fileName,
                           bool per_material=false,
                           int options = 0) const;
  void ExportCellsToVTK(const char* baseName) const;

  //02
  void BuildFaceHistogramInfo(double master_tolerance=100.0, double slave_tolerance=1.1);
  size_t NumberOfFaceHistogramBins();
  size_t MapFaceHistogramBins(size_t num_face_dofs);
  size_t GetFaceHistogramBinDOFSize(size_t category);
  bool IsCellLocal(uint64_t cell_global_index) const;
  bool IsCellBndry(uint64_t cell_global_index) const;

  void FindAssociatedVertices(chi_mesh::CellFace& cur_face,
                              std::vector<short>& dof_mapping) const;

  chi_mesh::Vector3
  ComputeCentroidFromListOfNodes(const std::vector<uint64_t>& list) const;

  std::vector<std::unique_ptr<chi_mesh::Cell>> GetGhostCells();

  ChiMPICommunicatorSet& GetCommunicator();

  size_t GetGlobalNumberOfCells() const;

  std::vector<uint64_t> GetDomainUniqueBoundaryIDs() const;

  size_t CountCellsInLogicalVolume(chi_mesh::LogicalVolume& log_vol) const;
};

#endif //CHI_MESHCONTINUUM_H_