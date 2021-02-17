#ifndef SPATIAL_DISCRETIZATION_FV_H
#define SPATIAL_DISCRETIZATION_FV_H

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "CellViews/fv_cellbase.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

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
  std::vector<chi_mesh::Cell*> neighbor_cells;
  std::vector<CellFVValues*> neighbor_cell_fv_views;

public:
  int              fv_local_block_address = 0;
  std::vector<int> locJ_block_address;
  std::vector<int> locJ_block_size;

private:
  explicit
  SpatialDiscretization_FV(chi_mesh::MeshContinuumPtr in_grid);

public:
  //prevent anything else other than a shared pointer
  static
  std::shared_ptr<SpatialDiscretization_FV>
  New(chi_mesh::MeshContinuumPtr in_grid)
  { return std::shared_ptr<SpatialDiscretization_FV>(
    new SpatialDiscretization_FV(in_grid));}

  //01
  void PreComputeCellSDValues(chi_mesh::MeshContinuumPtr grid) override;
  void PreComputeNeighborCellSDValues(chi_mesh::MeshContinuumPtr grid);

  CellFVValues* MapFeView(int cell_local_index);
  CellFVValues* MapNeighborFeView(int cell_global_index);

  //02 node ordering
  void OrderNodes(chi_mesh::MeshContinuumPtr grid);

  //03 sparsity
  void BuildSparsityPattern(chi_mesh::MeshContinuumPtr grid,
                            std::vector<int>& nodal_nnz_in_diag,
                            std::vector<int>& nodal_nnz_off_diag,
                            chi_math::UnknownManager& unknown_manager) override;

  //04a mappings
  int MapDOF(chi_mesh::Cell& cell);
  int MapDOFLocal(chi_mesh::Cell& cell);

  int MapDOF(chi_mesh::Cell& cell,
             chi_math::UnknownManager& unknown_manager,
             unsigned int unknown_id,
             unsigned int component=0);

  int MapDOFLocal(chi_mesh::Cell& cell,
                  chi_math::UnknownManager& unknown_manager,
                  unsigned int unknown_id,
                  unsigned int component=0);


  //04b utils
  size_t GetNumLocalDOFs(chi_mesh::MeshContinuumPtr grid,
                         chi_math::UnknownManager& unknown_manager) override;
  size_t GetNumGlobalDOFs(chi_mesh::MeshContinuumPtr grid,
                          chi_math::UnknownManager& unknown_manager) override;
  size_t GetNumGhostDOFs(chi_mesh::MeshContinuumPtr grid,
                         chi_math::UnknownManager& unknown_manager);

  std::vector<int> GetGhostDOFIndices(chi_mesh::MeshContinuumPtr grid,
                                      chi_math::UnknownManager& unknown_manager,
                                      unsigned int unknown_id=0);
  void LocalizePETScVector(Vec petsc_vector,
                           std::vector<double>& local_vector,
                           chi_math::UnknownManager& unknown_manager)
                           override;
};


#endif