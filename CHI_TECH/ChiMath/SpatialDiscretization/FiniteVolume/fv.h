#ifndef _chi_discretization_fv_h
#define _chi_discretization_fv_h

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "CellViews/fv_cellbase.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

//###################################################################
/**Spatial discretizations supporting Finite Volume representations.
 * */
class SpatialDiscretization_FV : public SpatialDiscretization
{
private:
  std::vector<CellFVView*> cell_fv_views;

private:
  bool               mapping_initialized;
  std::vector<bool>  cell_view_added_flags;

private:
  std::vector<chi_mesh::Cell*> neighbor_cells;
  std::vector<CellFVView*> neighbor_cell_fv_views;

public:
  SpatialDiscretization_FV(int dim=0,
                           chi_math::SpatialDiscretizationType sd_method =
                           chi_math::SpatialDiscretizationType::FINITE_VOLUME);

  //01
  void AddViewOfLocalContinuum(chi_mesh::MeshContinuum* grid) override;
  void AddViewOfNeighborContinuums(chi_mesh::MeshContinuum* grid);

  CellFVView* MapFeView(int cell_local_index);
  CellFVView* MapNeighborFeView(int cell_global_index);

  //02 node ordering
  void ReOrderNodes(chi_mesh::MeshContinuum* grid);

  //03 sparsity
  void BuildSparsityPattern(chi_mesh::MeshContinuum* grid,
                            std::vector<int>& nodal_nnz_in_diag,
                            std::vector<int>& nodal_nnz_off_diag,
                            chi_math::UnknownManager* unknown_manager=nullptr);

  //04a mappings
  int MapDOF(chi_mesh::Cell* cell);
  int MapDOFLocal(chi_mesh::Cell* cell);

  int MapDOF(chi_mesh::Cell* cell,
             chi_math::UnknownManager* unknown_manager,
             unsigned int unknown_id,
             unsigned int component=0);

  int MapDOFLocal(chi_mesh::Cell* cell,
                  chi_math::UnknownManager* unknown_manager,
                  unsigned int unknown_id,
                  unsigned int component=0);


  //04b utils
  unsigned int GetNumLocalDOFs(chi_mesh::MeshContinuum* grid,
                               chi_math::UnknownManager* unknown_manager=nullptr);
  unsigned int GetNumGlobalDOFs(chi_mesh::MeshContinuum* grid,
                                chi_math::UnknownManager* unknown_manager=nullptr);
  unsigned int GetNumGhostDOFs(chi_mesh::MeshContinuum* grid,
                               chi_math::UnknownManager* unknown_manager);

  std::vector<int> GetGhostDOFIndices(chi_mesh::MeshContinuum* grid,
                                      chi_math::UnknownManager* unknown_manager,
                                      unsigned int unknown_id=0);
  void LocalizePETScVector(Vec petsc_vector,
                           std::vector<double>& local_vector,
                           chi_math::UnknownManager* unknown_manager)
                           override;
};


#endif