#include "FLUDSCommonData.h"

namespace chi_mesh::sweep_management
{

FLUDSCommonData::FLUDSCommonData(
  const SPDS& spds,
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
  : spds_(spds), grid_nodal_mappings_(grid_nodal_mappings)
{
}

const SPDS& FLUDSCommonData::GetSPDS() const
{
  return spds_;
}

} // namespace chi_mesh::sweep_management