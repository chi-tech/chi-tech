#ifndef CHI_MESHCONTINUUM_H_
#define CHI_MESHCONTINUUM_H_

#include <memory>
#include <array>

#include "../chi_mesh.h"
#include "chi_meshcontinuum_localcellhandler.h"
#include "chi_meshcontinuum_globalcellhandler.h"
#include "chi_meshcontinuum_vertexhandler.h"

#include "chi_mpi.h"

namespace chi_data_types
{
  template<typename T>
  class NDArray;
}//namespace chi_data_types

//######################################################### Class Definition
/**Stores the relevant information for completely defining a computational
 * domain. */
class chi_mesh::MeshContinuum
{
private:
  std::vector<std::unique_ptr<chi_mesh::Cell>> native_cells;  ///< Actual native cells
  std::vector<std::unique_ptr<chi_mesh::Cell>> foreign_cells; ///< Locally stored ghosts

  std::map<uint64_t,uint64_t> global_cell_id_to_native_id_map;
  std::map<uint64_t,uint64_t> global_cell_id_to_foreign_id_map;

  uint64_t global_vertex_count=0;

public:
  VertexHandler                  vertices;
  LocalCellHandler               local_cells;
  GlobalCellHandler              cells;
  std::vector<uint64_t>          local_cell_glob_indices;

private:
  bool                           face_histogram_available = false;
  bool                           communicators_available  = false;

  //Pair.first is the max dofs-per-face for the category and Pair.second
  //is the number of faces in this category
  std::vector<std::pair<size_t,size_t>> face_categories;

  chi_objects::ChiMPICommunicatorSet communicator_set;

  MeshAttributes attributes = NONE;

  struct
  {
    size_t Nx = 0;
    size_t Ny = 0;
    size_t Nz = 0;
  }ortho_attributes;

public:
  MeshContinuum() :
    local_cells(native_cells),
    cells(local_cell_glob_indices,
          native_cells,
          foreign_cells,
          global_cell_id_to_native_id_map,
          global_cell_id_to_foreign_id_map)
  {
  }

  void SetGlobalVertexCount(const uint64_t count) {global_vertex_count = count;}
  uint64_t GetGlobalVertexCount() const {return global_vertex_count;}

  static
  std::shared_ptr<MeshContinuum> New()
  { return std::make_shared<MeshContinuum>();}

  /**Method to be called if cells and nodes have been transferred
   * to another grid.*/
  void ClearCellReferences()
  {
    native_cells.clear();
    foreign_cells.clear();
    global_cell_id_to_native_id_map.clear();
    global_cell_id_to_foreign_id_map.clear();
  }

  void ExportCellsToObj(const char* fileName,
                           bool per_material=false,
                           int options = 0) const;
  void ExportCellsToVTK(const char* baseName) const;

  void BuildFaceHistogramInfo(double master_tolerance=100.0, double slave_tolerance=1.1);
  size_t NumberOfFaceHistogramBins();
  size_t MapFaceHistogramBins(size_t num_face_dofs);
  size_t GetFaceHistogramBinDOFSize(size_t category);
  bool IsCellLocal(uint64_t cell_global_index) const;
  bool IsCellBndry(uint64_t cell_global_index) const;

  void FindAssociatedVertices(const chi_mesh::CellFace& cur_face,
                              std::vector<short>& dof_mapping) const;
  static size_t MapCellFace(const chi_mesh::Cell& cur_cell,
                     const chi_mesh::Cell& adj_cell,
                     unsigned int f) ;

  chi_mesh::Vector3
  ComputeCentroidFromListOfNodes(const std::vector<uint64_t>& list) const;

  chi_objects::ChiMPICommunicatorSet& GetCommunicator();

  size_t GetGlobalNumberOfCells() const;

  std::vector<uint64_t> GetDomainUniqueBoundaryIDs() const;

  size_t CountCellsInLogicalVolume(const chi_mesh::LogicalVolume& log_vol) const;

  MeshAttributes Attributes() const {return attributes;}

  std::array<size_t,3> GetIJKInfo() const;
  chi_data_types::NDArray<uint64_t> MakeIJKToGlobalIDMapping() const;
  std::vector<chi_mesh::Vector3> MakeCellOrthoSizes() const;

private:
  friend class chi_mesh::VolumeMesher;
  void SetAttributes(MeshAttributes new_attribs,
                     std::array<size_t,3> ortho_Nis={0,0,0})
  {
    attributes = attributes | new_attribs;
    ortho_attributes = {ortho_Nis[0],ortho_Nis[1],ortho_Nis[2]};
  }
};

#endif //CHI_MESHCONTINUUM_H_