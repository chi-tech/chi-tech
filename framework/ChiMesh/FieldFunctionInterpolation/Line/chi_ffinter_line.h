#ifndef CHI_FFINTER_LINE_H
#define CHI_FFINTER_LINE_H

#include "../chi_ffinterpolation.h"

#include <petscksp.h>

namespace chi_mesh
{
  struct FieldFunctionContext
  {
    std::shared_ptr<chi_physics::FieldFunction>    ref_ff;
    std::vector<double>            interpolation_points_values;
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
  int                number_of_points = 2;
  chi_mesh::Vector3  pi, pf;

  std::vector<std::vector<double>>   custom_arrays;
  std::vector<chi_mesh::Vector3>     interpolation_points;
  std::vector<FieldFunctionContext>  ff_contexts;

private:
  double                         delta_d = 1.0;

public:
  FieldFunctionInterpolationLine() :
    FieldFunctionInterpolation(ff_interpolation::Type::LINE)
  {  }
  //01
  void Initialize() override;

  //02
  void Execute() override;

public:
  std::string GetDefaultFileBaseName() const override
  {return "ZLFFI";}
  void ExportPython(std::string base_name) override;
};

#endif