#include "BooleanLogicalVolume.h"

#include "ChiObject/object_maker.h"

namespace chi_mesh
{

RegisterChiObject(chi_mesh, BooleanLogicalVolume);

chi_objects::InputParameters BooleanLogicalVolume::GetInputParameters()
{
  chi_objects::InputParameters params = LogicalVolume::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription(
  "\\defgroup chi_mesh__BooleanLogicalVolume BooleanLogicalVolume\n"
  "\\ingroup LuaLogicVolumes");
  // clang-format on

  params.AddRequiredParameterArray("parts", "Array of combinatorial logic");

  return params;
}

BooleanLogicalVolume::BooleanLogicalVolume(
  const chi_objects::InputParameters& params)
  : LogicalVolume(params)
{
  const auto& input_parts = params.GetParam("parts");
  input_parts.RequireBlockTypeIs(chi_objects::ParameterBlockType::ARRAY);

  for (size_t p = 0; p < input_parts.NumParameters(); ++p)
  {
    const auto& part = input_parts.GetParam(p);
    part.RequireBlockTypeIs(chi_objects::ParameterBlockType::BLOCK);

    part.RequireParameter("op");
    part.RequireParameter("lv");

    const size_t lv_handle = part.GetParamValue<size_t>("lv");
    auto lv_ptr = chi::GetStackItemPtrAsType<LogicalVolume>(
      chi::object_stack, lv_handle, __FUNCTION__);

    parts.emplace_back(part.GetParamValue<bool>("op"), lv_ptr);
  }
}

bool BooleanLogicalVolume::Inside(const chi_mesh::Vector3& point) const
{
  for (const auto& part : parts)
  {
    if (not(part.first && part.second->Inside(point))) return false;
  }

  return true;
}

} // namespace chi_mesh