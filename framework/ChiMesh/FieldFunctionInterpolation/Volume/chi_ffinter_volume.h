#ifndef CHI_FFINTER_VOLUME_H
#define CHI_FFINTER_VOLUME_H

#include "../chi_ffinterpolation.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#include <petscksp.h>

//###################################################################
/**Volume-wise field function interpolation.
 *
 * This interpolator allows the user to obtain quantities by logical
 * volume. If no logical volume is assigned to the method it will
 * default to operating over the entire volume.\n
 * \n
 * The method also supports a few primitive operations:
 *  - OP_VOLUME_AVG. Obtains the volume average of the field function
 *    of interest.
 *  - OP_VOLUME_SUM. Obtains the volume integral of the field function
 *    of interest.*/
class chi_mesh::FieldFunctionInterpolationVolume :
  public chi_mesh::FieldFunctionInterpolation
{
public:
  std::shared_ptr<chi_mesh::LogicalVolume> logical_volume = nullptr;
  ff_interpolation::Operation op_type = ff_interpolation::Operation::OP_SUM;
  std::string op_lua_func;
  double op_value = 0.0;

private:
  std::vector<uint64_t> cell_local_ids_inside_logvol;

public:
  FieldFunctionInterpolationVolume() :
    FieldFunctionInterpolation(ff_interpolation::Type::VOLUME)
  {  }

  //01
  void Initialize() override;

  //02
  void Execute() override;

  void CFEMInterpolate(Vec field, std::vector<uint64_t> &mapping);
  void PWLDInterpolate(std::vector<double>& field, std::vector<uint64_t> &mapping);

  double CallLuaFunction(double ff_value, int mat_id);

  std::string GetDefaultFileBaseName() const override
  {return "ZVFFI";}
  void ExportPython(std::string base_name) override {}
};


#endif
