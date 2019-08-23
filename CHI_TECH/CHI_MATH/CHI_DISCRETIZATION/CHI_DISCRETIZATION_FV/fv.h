#ifndef _chi_discretization_fv_h
#define _chi_discretization_fv_h

#include "../chi_discretization.h"
#include "CellViews/fv_cellbase.h"

//###################################################################
/**Spatial discretizations supporting Finite Volume representations.
 * */
class CHI_DISCRETIZATION_FV : public CHI_DISCRETIZATION
{
private:
  std::vector<CellFVView*> cell_fv_views;
  std::vector<int>         cell_fv_views_mapping;

private:
  bool mapping_initialized;


public:
  CHI_DISCRETIZATION_FV(int dim=0);

  void AddViewOfLocalContinuum(
    chi_mesh::MeshContinuum* vol_continuum,
    int num_cells,
    int* cell_indices);

  CellFVView* MapFeView(int cell_glob_index);
};


#endif