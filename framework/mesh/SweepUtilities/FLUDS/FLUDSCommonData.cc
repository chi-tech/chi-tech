#include "FLUDSCommonData.h"

namespace chi_mesh::sweep_management
{

FLUDSCommonData::FLUDSCommonData(
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
  : grid_nodal_mappings_(grid_nodal_mappings)
{
}

} // namespace chi_mesh::sweep_management