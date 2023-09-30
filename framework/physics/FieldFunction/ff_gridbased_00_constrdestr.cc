#include "fieldfunction_gridbased.h"

#include "math/SpatialDiscretization/FiniteVolume/fv.h"
#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"
#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "mesh/MeshHandler/chi_meshhandler.h"

#include "ChiObjectFactory.h"

namespace chi_physics
{

RegisterChiObject(chi_physics, FieldFunctionGridBased);

chi::InputParameters FieldFunctionGridBased::GetInputParameters()
{
  chi::InputParameters params = FieldFunction::GetInputParameters();

  params.SetDocGroup("DocFieldFunction");

  params.AddRequiredParameter<std::string>(
    "sdm_type", "The spatial discretization type to be used");

  params.AddOptionalParameter(
    "initial_value", 0.0, "The initial value to assign to the field function");

  params.AddOptionalParameter(
    "quadrature_order",
    0,
    "If supplied, will overwrite the default for the "
    "specific discretization-coordinate system combination.");

  params.AddOptionalParameter("coordinate_system",
                              "cartesian",
                              "Coordinate system to apply to element mappings");

  using namespace chi_data_types;
  params.ConstrainParameterRange(
    "sdm_type", AllowableRangeList::New({"FV", "PWLC", "PWLD"}));
  params.ConstrainParameterRange(
    "coordinate_system",
    AllowableRangeList::New({"cartesian", "rz", "1d_spherical"}));

  return params;
}

// ##################################################################
/**ObjectMaker based constructor.*/
FieldFunctionGridBased::FieldFunctionGridBased(
  const chi::InputParameters& params)
  : FieldFunction(params),
    sdm_(MakeSpatialDiscretization(params)),
    ghosted_field_vector_(MakeFieldVector(*sdm_, GetUnknownManager())),
    local_grid_bounding_box_(
      chi_mesh::GetCurrentHandler().GetGrid()->GetLocalBoundingBox())
{
  ghosted_field_vector_->Set(params.GetParamValue<double>("initial_value"));
}

// ##################################################################
FieldFunctionGridBased::FieldFunctionGridBased(
  const std::string& text_name,
  chi_math::SDMPtr& discretization_ptr,
  chi_math::Unknown unknown)
  : FieldFunction(text_name, std::move(unknown)),
    sdm_(discretization_ptr),
    ghosted_field_vector_(MakeFieldVector(*sdm_, GetUnknownManager())),
    local_grid_bounding_box_(sdm_->Grid().GetLocalBoundingBox())
{
}

// ##################################################################
FieldFunctionGridBased::FieldFunctionGridBased(
  const std::string& text_name,
  chi_math::SDMPtr& sdm_ptr,
  chi_math::Unknown unknown,
  const std::vector<double>& field_vector)
  : FieldFunction(text_name, std::move(unknown)),
    sdm_(sdm_ptr),
    ghosted_field_vector_(MakeFieldVector(*sdm_, GetUnknownManager())),
    local_grid_bounding_box_(sdm_->Grid().GetLocalBoundingBox())
{
  ChiInvalidArgumentIf(
    field_vector.size() != ghosted_field_vector_->LocalSize(),
    "Constructor called with incompatible size field vector.");

  ghosted_field_vector_->Set(field_vector);
}

// ##################################################################
FieldFunctionGridBased::FieldFunctionGridBased(const std::string& text_name,
                                               chi_math::SDMPtr& sdm_ptr,
                                               chi_math::Unknown unknown,
                                               double field_value)
  : FieldFunction(text_name, std::move(unknown)),
    sdm_(sdm_ptr),
    ghosted_field_vector_(MakeFieldVector(*sdm_, GetUnknownManager())),
    local_grid_bounding_box_(sdm_->Grid().GetLocalBoundingBox())
{
  ghosted_field_vector_->Set(field_value);
}

// ##################################################################
/**Returns the spatial discretization method.*/
const chi_math::SpatialDiscretization&
FieldFunctionGridBased::GetSpatialDiscretization() const
{
  return *sdm_;
}

// ##################################################################
/**Returns a read-only reference to the locally stored field data.*/
const std::vector<double>& FieldFunctionGridBased::FieldVectorRead() const
{
  return ghosted_field_vector_->LocalSTLData();
}
/**Returns a reference to the locally stored field data.*/
std::vector<double>& FieldFunctionGridBased::FieldVector()
{
  return ghosted_field_vector_->LocalSTLData();
}

// ##################################################################
/**Private method for creating the spatial discretization method.*/
chi_math::SDMPtr FieldFunctionGridBased::MakeSpatialDiscretization(
  const chi::InputParameters& params)
{
  const auto& user_params = params.ParametersAtAssignment();
  const auto& grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto sdm_type = params.GetParamValue<std::string>("sdm_type");

  typedef chi_math::SpatialDiscretization_FV FV;
  typedef chi_math::SpatialDiscretization_PWLC PWLC;
  typedef chi_math::SpatialDiscretization_PWLD PWLD;

  if (sdm_type == "FiniteVolume") return FV::New(*grid_ptr);

  chi_math::finite_element::SetupFlags setup_flags =
    chi_math::finite_element::NO_FLAGS_SET;
  chi_math::QuadratureOrder q_order = chi_math::QuadratureOrder::SECOND;
  chi_math::CoordinateSystemType coordinate_system_type =
    chi_math::CoordinateSystemType::CARTESIAN;

  if (user_params.Has("quadrature_order"))
  {
    const uint32_t max_order =
      static_cast<uint32_t>(chi_math::QuadratureOrder::FORTYTHIRD);
    const uint32_t q_order_int =
      params.GetParamValue<uint32_t>("quadrature_order");
    ChiInvalidArgumentIf(q_order_int > max_order,
                         "Invalid quadature point order " +
                           std::to_string(q_order_int));
    q_order = static_cast<chi_math::QuadratureOrder>(q_order_int);
  }

  if (user_params.Has("coordinate_system"))
  {
    const auto cs = params.GetParamValue<std::string>("coordinate_system");

    if (cs == "cartesian")
      coordinate_system_type = chi_math::CoordinateSystemType::CARTESIAN;
    if (cs == "rz")
      coordinate_system_type = chi_math::CoordinateSystemType::CYLINDRICAL;
    if (cs == "1d_spherical")
      coordinate_system_type = chi_math::CoordinateSystemType::SPHERICAL;
  }

  if (sdm_type == "PWLC")
    return PWLC::New(*grid_ptr, setup_flags, q_order, coordinate_system_type);
  if (sdm_type == "PWLD")
    return PWLD::New(*grid_ptr, setup_flags, q_order, coordinate_system_type);

  // If not returned by now
  ChiInvalidArgument("Unsupported sdm_type \"" + sdm_type + "\"");
}

/**Private method for creating the field vector.*/
std::unique_ptr<chi_math::GhostedParallelSTLVector>
FieldFunctionGridBased::MakeFieldVector(
  const chi_math::SpatialDiscretization& discretization,
  const chi_math::UnknownManager& uk_man)
{
  auto field = std::make_unique<chi_math::GhostedParallelSTLVector>(
    discretization.GetNumLocalDOFs(uk_man),
    discretization.GetNumGlobalDOFs(uk_man),
    discretization.GetGhostDOFIndices(uk_man),
    Chi::mpi.comm);

  return field;
}

} // namespace chi_physics