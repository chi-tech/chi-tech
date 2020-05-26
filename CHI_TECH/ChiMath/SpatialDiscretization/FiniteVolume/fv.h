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


public:
  SpatialDiscretization_FV(int dim=0);

  void AddViewOfLocalContinuum(chi_mesh::MeshContinuum* grid) override;

  int MapDOF(chi_mesh::Cell* cell,
             chi_math::UnknownManager* unknown_manager,
             unsigned int unknown_id,
             unsigned int component=0);
  int MapDOF(int cell_global_id,
             chi_math::UnknownManager* unknown_manager,
             unsigned int unknown_id,
             unsigned int component=0);

  int MapDOFLocal(chi_mesh::Cell* cell,
                  chi_math::UnknownManager* unknown_manager,
                  unsigned int unknown_id,
                  unsigned int component=0);

  CellFVView* MapFeView(int cell_local_index);

  void BuildSparsityPattern(chi_mesh::MeshContinuum* grid,
                            std::vector<int>& nodal_nnz_in_diag,
                            std::vector<int>& nodal_nnz_off_diag,
                            chi_math::UnknownManager* unknown_manager=nullptr);

  unsigned int GetNumLocalDOFs(chi_mesh::MeshContinuum* grid,
                               chi_math::UnknownManager* unknown_manager);
  unsigned int GetNumGlobalDOFs(chi_mesh::MeshContinuum* grid,
                                chi_math::UnknownManager* unknown_manager);
  unsigned int GetNumGhostDOFs(chi_mesh::MeshContinuum* grid,
                               chi_math::UnknownManager* unknown_manager);

  std::vector<int> GetGhostDOFIndices(chi_mesh::MeshContinuum* grid,
                                      chi_math::UnknownManager* unknown_manager,
                                      unsigned int unknown_id=0);
};


#endif