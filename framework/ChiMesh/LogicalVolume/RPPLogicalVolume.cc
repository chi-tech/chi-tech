#include "RPPLogicalVolume.h"

#include "ChiObject/object_maker.h"

namespace chi_mesh
{

RegisterChiObject(chi_mesh, RPPLogicalVolume);

chi_objects::InputParameters RPPLogicalVolume::GetInputParameters()
{
  chi_objects::InputParameters params = LogicalVolume::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription(
  "\\defgroup chi_mesh__RPPLogicalVolume RPPLogicalVolume\n"
  "\\ingroup LuaLogicVolumes");
  // clang-format on

  params.AddOptionalParameter("xmin", 0.0, "X-min of the volume");
  params.AddOptionalParameter("xmax", 1.0, "X-max of the volume");
  params.AddOptionalParameter("ymin", 0.0, "Y-min of the volume");
  params.AddOptionalParameter("ymax", 1.0, "Y-max of the volume");
  params.AddOptionalParameter("zmin", 0.0, "Z-min of the volume");
  params.AddOptionalParameter("zmax", 1.0, "Z-max of the volume");

  params.AddOptionalParameter(
    "infx", false, "Flag, when true, will ignore xmin and xmax.");
  params.AddOptionalParameter(
    "infy", false, "Flag, when true, will ignore ymin and ymax.");
  params.AddOptionalParameter(
    "infz", false, "Flag, when true, will ignore zmin and zmax.");

  return params;
}

RPPLogicalVolume::RPPLogicalVolume(const chi_objects::InputParameters& params)
  : LogicalVolume(params),
    xmin_(params.GetParamValue<double>("xmin")),
    xmax_(params.GetParamValue<double>("xmax")),
    ymin_(params.GetParamValue<double>("ymin")),
    ymax_(params.GetParamValue<double>("ymax")),
    zmin_(params.GetParamValue<double>("zmin")),
    zmax_(params.GetParamValue<double>("zmax"))
{
}

bool RPPLogicalVolume::Inside(const chi_mesh::Vector3& point) const
{
  if ((point.x <= xmax_) && (point.x >= xmin_) && (point.y <= ymax_) &&
      (point.y >= ymin_) && (point.z <= zmax_) && (point.z >= zmin_))
  {
    return true;
  }
  else
    return false;
}

} // namespace chi_mesh