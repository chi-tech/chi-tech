#include "BooleanLogicalVolume.h"

#include "ChiObject/object_maker.h"

namespace chi_mesh
{

chi_objects::InputParameters BooleanLogicalVolumeArgumentPair();

RegisterChiObject(chi_mesh, BooleanLogicalVolume);
RegisterSyntaxBlock(chi_mesh,
                    BooleanLogicalVolumeArgumentPair,
                    BooleanLogicalVolumeArgumentPair);

chi_objects::InputParameters BooleanLogicalVolume::GetInputParameters()
{
  chi_objects::InputParameters params = LogicalVolume::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription(
  "\\defgroup chi_mesh__BooleanLogicalVolume BooleanLogicalVolume\n"
  "\\ingroup LuaLogicVolumes\n");
  // clang-format on

  params.AddRequiredParameterArray(
    "parts",
    "Array of combinatorial logic each entry has the following required params "
    "<TT>chi_mesh::BooleanLogicalVolumeArgumentPair</TT>"
    "$(chi_mesh::BooleanLogicalVolumeArgumentPair$)");

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

    auto part_params = BooleanLogicalVolumeArgumentPair();

    part_params.AssignParameters(part);

    const size_t lv_handle = part_params.GetParamValue<size_t>("lv");
    auto lv_ptr = chi::GetStackItemPtrAsType<LogicalVolume>(
      chi::object_stack, lv_handle, __FUNCTION__);

    parts.emplace_back(part_params.GetParamValue<bool>("op"), lv_ptr);
  }
}

chi_objects::InputParameters BooleanLogicalVolumeArgumentPair()
{
  chi_objects::InputParameters params;

  // clang-format off
  params.SetGeneralDescription(
  "\\defgroup chi_mesh__BooleanLogicalVolumeArgumentPair chi_mesh.BooleanLogicalVolumeArgumentPair\n"
  "\\ingroup chi_mesh__BooleanLogicalVolume");
  // clang-format on

  params.AddRequiredParameter<bool>(
    "op",
    "Boolean value indicating the volume sense. True means inside, False means "
    "outside");
  params.AddRequiredParameter<size_t>("lv", "Handle to a logical volume.");

  return params;
}

bool BooleanLogicalVolume::Inside(const chi_mesh::Vector3& point) const
{
  for (const auto& part : parts)
  {
    if (part.first != part.second->Inside(point)) return false;
  }

  return true;
}

} // namespace chi_mesh