#ifndef _chi_ffinter_line_h
#define _chi_ffinter_line_h

#include "../chi_ffinterpolation.h"
#include "../../chi_mesh.h"

#include <petscksp.h>

#define FFI_LINE_FIRSTPOINT  11
#define FFI_LINE_SECONDPOINT 12
#define FFI_LINE_NUMBEROFPOINTS 13

struct FieldFunctionContext
{
  chi_physics::FieldFunction*    ref_ff;
  std::vector<double>            interpolation_points_values;
  std::vector<int>               cfem_local_nodes_needed_unmapped;
  std::vector<int>               pwld_local_nodes_needed_unmapped;
  std::vector<int>               pwld_local_cells_needed_unmapped;
  std::vector<int>               interpolation_points_ass_cell;
};

//###################################################################
/** A line based interpolation function.*/
class chi_mesh::FieldFunctionInterpolationLine :
  public FieldFunctionInterpolation
{
public:
  int               number_of_points;
  chi_mesh::Vector  pi,pf;

private:
  double                         delta_d;
  std::vector<chi_mesh::Vector>  interpolation_points;
  std::vector<double>            interpolation_points_values;
  std::vector<int>               cfem_local_nodes_needed_unmapped;
  std::vector<int>               pwld_local_nodes_needed_unmapped;
  std::vector<int>               pwld_local_cells_needed_unmapped;
  std::vector<int>               interpolation_points_ass_cell;

  std::vector<FieldFunctionContext*> ff_contexts;
public:
  FieldFunctionInterpolationLine()
  {
    number_of_points = 2;
    delta_d = 1.0;
  }
  //01
  void Initialize();

  //02
  void Execute();
private:
  void CFEMInterpolate(Vec field, std::vector<int> &mapping);
  void PWLDInterpolate(std::vector<double>& field, std::vector<int> &mapping);

  void CFEMInterpolate(Vec field,
                       std::vector<int> &mapping,
                       FieldFunctionContext* ff_ctx);
  void PWLDInterpolate(std::vector<int> &mapping,
                       FieldFunctionContext* ff_ctx);
  void FVInterpolate(std::vector<int> &mapping,
                     FieldFunctionContext* ff_ctx);
public:
  void ExportPython(std::string base_name);
};

#endif