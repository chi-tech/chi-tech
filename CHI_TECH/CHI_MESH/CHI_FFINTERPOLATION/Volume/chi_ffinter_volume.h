#ifndef _chi_ffinter_volume_h
#define _chi_ffinter_volume_h

#include "../chi_ffinterpolation.h"
#include <CHI_MESH/CHI_LOGICALVOLUME/chi_mesh_logicalvolume.h>

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
  chi_mesh::LogicalVolume* logical_volume;
  int op_type;
  double op_value;

private:
  std::vector<int>                    cfem_local_nodes_needed_unmapped;
  std::vector<int>                    pwld_local_nodes_needed_unmapped;
  std::vector<int>                    pwld_local_cells_needed_unmapped;

public:
  FieldFunctionInterpolationVolume()
  {
    logical_volume = nullptr;
    op_type = OP_SUM;
    op_value = 0.0;
  }

  //01
  void Initialize();

  //02
  void Execute();

  void CFEMInterpolate(Vec field, std::vector<int> &mapping);
  void PWLDInterpolate(std::vector<double>& field, std::vector<int> &mapping);

  void ExportPython(std::string base_name)
  {

  }

};


#endif
