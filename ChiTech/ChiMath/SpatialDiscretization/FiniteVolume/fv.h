#ifndef SPATIAL_DISCRETIZATION_FV_H
#define SPATIAL_DISCRETIZATION_FV_H

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "CellViews/fv_cellbase.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <map>

//###################################################################
/**Spatial discretizations supporting Finite Volume representations.
 * */
class SpatialDiscretization_FV : public SpatialDiscretization
{
private:
  std::vector<CellFVValues*> cell_fv_views;

private:
  bool               mapping_initialized;
  std::vector<bool>  cell_view_added_flags;

private:
  std::map<uint64_t, std::unique_ptr<chi_mesh::Cell>> neighbor_cells;
  std::map<uint64_t, CellFVValues*> neighbor_cell_fv_views;

public:
  int              fv_local_block_address = 0;
  std::vector<int> locJ_block_address;
  std::vector<int> locJ_block_size;

private:
  explicit
  SpatialDiscretization_FV(chi_mesh::MeshContinuumPtr& in_grid,
                           chi_math::CoordinateSystemType in_cs_type);

public:
  virtual ~SpatialDiscretization_FV() = default;
  //prevent anything else other than a shared pointer
  static
  std::shared_ptr<SpatialDiscretization_FV>
  New(chi_mesh::MeshContinuumPtr& in_grid,
      chi_math::CoordinateSystemType in_cs_type =
      chi_math::CoordinateSystemType::CARTESIAN)
  { return std::shared_ptr<SpatialDiscretization_FV>(
    new SpatialDiscretization_FV(in_grid, in_cs_type));}

  //01
  void PreComputeCellSDValues() override;
  void PreComputeNeighborCellSDValues();

  CellFVValues* MapFeView(uint64_t cell_local_index);
  CellFVValues* MapNeighborFeView(uint64_t cell_global_index);

  //02 node ordering
  void OrderNodes();

  //03 sparsity
  void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                            std::vector<int64_t>& nodal_nnz_off_diag,
                            chi_math::UnknownManager& unknown_manager) override;

  //04a mappings
  int64_t MapDOF(const chi_mesh::Cell& cell,
                 unsigned int node,
                 const chi_math::UnknownManager& unknown_manager,
                 unsigned int unknown_id,
                 unsigned int component) const override;

  int64_t MapDOFLocal(const chi_mesh::Cell& cell,
                  unsigned int node,
                  const chi_math::UnknownManager& unknown_manager,
                  unsigned int unknown_id,
                  unsigned int component) const override;

  int64_t MapDOF(const chi_mesh::Cell& cell, unsigned int node) const override;
  int64_t MapDOFLocal(const chi_mesh::Cell& cell, unsigned int node) const override;

  //04b utils
  size_t GetNumLocalDOFs(chi_math::UnknownManager& unknown_manager) override;
  size_t GetNumGlobalDOFs(chi_math::UnknownManager& unknown_manager) override;
  size_t GetNumGhostDOFs(chi_mesh::MeshContinuumPtr grid,
                         chi_math::UnknownManager& unknown_manager);

  size_t GetCellNumNodes(const chi_mesh::Cell& cell) const override
  {return 1;}

  std::vector<chi_mesh::Vector3>
    GetCellNodeLocations(const chi_mesh::Cell& cell) const override
  {
    std::vector<chi_mesh::Vector3> node_locations(1,cell.centroid);

    return node_locations;
  }

  std::vector<int> GetGhostDOFIndices(chi_mesh::MeshContinuumPtr grid,
                                      chi_math::UnknownManager& unknown_manager,
                                      unsigned int unknown_id=0);
  void LocalizePETScVector(Vec petsc_vector,
                           std::vector<double>& local_vector,
                           chi_math::UnknownManager& unknown_manager)
                           override;
};


#endif