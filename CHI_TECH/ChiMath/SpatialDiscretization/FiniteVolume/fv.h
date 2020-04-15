#ifndef _chi_discretization_fv_h
#define _chi_discretization_fv_h

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "CellViews/fv_cellbase.h"

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

  CellFVView* MapFeView(int cell_local_index);

  void BuildSparsityPattern(chi_mesh::MeshContinuum* grid,
                            std::vector<int>& nodal_nnz_in_diag,
                            std::vector<int>& nodal_nnz_off_diag);
};


#endif