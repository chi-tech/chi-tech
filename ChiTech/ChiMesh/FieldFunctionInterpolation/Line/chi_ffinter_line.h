#ifndef CHI_FFINTER_LINE_H
#define CHI_FFINTER_LINE_H

#include "../chi_ffinterpolation.h"
#include "../../chi_mesh.h"

#include <petscksp.h>

#define FFI_LINE_FIRSTPOINT  11
#define FFI_LINE_SECONDPOINT 12
#define FFI_LINE_NUMBEROFPOINTS 13
#define FFI_LINE_CUSTOM_ARRAY 14

namespace chi_mesh
{
  struct FieldFunctionContext
  {
    std::shared_ptr<chi_physics::FieldFunction>    ref_ff;
    std::vector<double>            interpolation_points_values;
    std::vector<int>               cfem_local_nodes_needed_unmapped;
    std::vector<int>               cfem_local_cells_needed_unmapped;
    std::vector<int>               pwld_local_nodes_needed_unmapped;
    std::vector<int>               pwld_local_cells_needed_unmapped;
    std::vector<uint64_t>          interpolation_points_ass_cell;
    std::vector<bool>              interpolation_points_has_ass_cell;
  };
}

//###################################################################
/** A line based interpolation function.*/
class chi_mesh::FieldFunctionInterpolationLine :
  public FieldFunctionInterpolation
{
public:
  int               number_of_points;
  chi_mesh::Vector3  pi,pf;

  std::vector<std::vector<double>> custom_arrays;
  std::vector<chi_mesh::Vector3>  interpolation_points;
  std::vector<FieldFunctionContext*> ff_contexts;

private:
  double                         delta_d;

public:
  FieldFunctionInterpolationLine()
  {
    number_of_points = 2;
    delta_d = 1.0;
  }
  //01
  void Initialize() override;

  //02
  void Execute() override;
private:
  void CFEMInterpolate(Vec field,
                       std::vector<uint64_t>& mapping,
                       FieldFunctionContext* ff_ctx);
  void PWLDInterpolate(std::vector<uint64_t>& mapping,
                       FieldFunctionContext* ff_ctx);
  void FVInterpolate(std::vector<uint64_t>& mapping,
                     FieldFunctionContext* ff_ctx);
public:
  std::string GetDefaultFileBaseName() const override
  {return "ZLFFI";}
  void ExportPython(std::string base_name) override;
};

#endif