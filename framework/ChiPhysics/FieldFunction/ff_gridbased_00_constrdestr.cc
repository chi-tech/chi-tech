#include "fieldfunction_gridbased.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "ChiObject/object_maker.h"

namespace chi_physics
{

RegisterChiObject(chi_physics, FieldFunctionGridBased);

chi::InputParameters FieldFunctionGridBased::GetInputParameters()
{
  chi::InputParameters params = FieldFunction::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription(
  "\\defgroup chi_physics__FieldFunctionGridBased "
  " chi_physics.FieldFunctionGridBased\n"
  "\\ingroup DocFieldFunction");
  // clang-format on

  params.AddRequiredParameter<std::string>(
    "sdm_type", "The spatial discretization type to be used");

  params.AddOptionalParameter(
    "initial_value", 0.0, "The initial value to assign to the field function");

  return params;
}

// ##################################################################
/**ObjectMaker based constructor.*/
FieldFunctionGridBased::FieldFunctionGridBased(
  const chi::InputParameters& params)
  : FieldFunction(params),
    local_grid_bounding_box_(
      chi_mesh::GetCurrentHandler().GetGrid()->GetLocalBoundingBox())
{
  const auto& grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto sdm_type = params.GetParamValue<std::string>("sdm_type");

  typedef chi_math::SpatialDiscretization_FV FV;
  typedef chi_math::SpatialDiscretization_PWLC PWLC;
  typedef chi_math::SpatialDiscretization_PWLD PWLD;

  if (sdm_type == "FiniteVolume") sdm_ = FV::New(*grid_ptr);
  if (sdm_type == "PWLC") sdm_ = PWLC::New(*grid_ptr);
  if (sdm_type == "PWLD") sdm_ = PWLD::New(*grid_ptr);

  const size_t num_local_dofs = sdm_->GetNumLocalDOFs(UnkManager());
  field_vector_.assign(num_local_dofs,
                       params.GetParamValue<double>("initial_value"));

  vector_ghost_communicator_ = MakeGhostCommunicator();
}

// ##################################################################
FieldFunctionGridBased::FieldFunctionGridBased(const std::string& text_name,
                                               chi_math::SMDPtr& sdm_ptr,
                                               chi_math::Unknown unknown)
  : FieldFunction(text_name, std::move(unknown)),
    sdm_(sdm_ptr),
    local_grid_bounding_box_(sdm_->Grid().GetLocalBoundingBox())
{
  const size_t num_local_dofs = sdm_->GetNumLocalDOFs(UnkManager());
  field_vector_.assign(num_local_dofs, 0.0);

  vector_ghost_communicator_ = MakeGhostCommunicator();
}

// ##################################################################
FieldFunctionGridBased::FieldFunctionGridBased(const std::string& text_name,
                                               chi_math::SMDPtr& sdm_ptr,
                                               chi_math::Unknown unknown,
                                               std::vector<double> field_vector)
  : FieldFunction(text_name, std::move(unknown)),
    sdm_(sdm_ptr),
    local_grid_bounding_box_(sdm_->Grid().GetLocalBoundingBox())
{
  const std::string fname = __FUNCTION__;
  const size_t num_local_dofs = sdm_->GetNumLocalDOFs(UnkManager());
  if (field_vector.size() != num_local_dofs)
    throw std::logic_error(fname +
                           ": Constructor initialized with incompatible "
                           "size field vector.");

  field_vector_ = std::move(field_vector);

  vector_ghost_communicator_ = MakeGhostCommunicator();
}

// ##################################################################
FieldFunctionGridBased::FieldFunctionGridBased(const std::string& text_name,
                                               chi_math::SMDPtr& sdm_ptr,
                                               chi_math::Unknown unknown,
                                               double field_value)
  : FieldFunction(text_name, std::move(unknown)),
    sdm_(sdm_ptr),
    local_grid_bounding_box_(sdm_->Grid().GetLocalBoundingBox())
{
  const size_t num_local_dofs = sdm_->GetNumLocalDOFs(UnkManager());
  field_vector_.assign(num_local_dofs, field_value);

  vector_ghost_communicator_ = MakeGhostCommunicator();
}

// ##################################################################
/**Private method for creating the ghost communicator.*/
std::shared_ptr<chi_math::VectorGhostCommunicator>
chi_physics::FieldFunctionGridBased::MakeGhostCommunicator()
{
  const size_t num_local_dofs = sdm_->GetNumLocalDOFs(UnkManager());
  const size_t num_globl_dofs = sdm_->GetNumGlobalDOFs(UnkManager());
  const std::vector<int64_t> ghost_ids = sdm_->GetGhostDOFIndices(UnkManager());

  return std::make_shared<chi_math::VectorGhostCommunicator>(
    num_local_dofs, num_globl_dofs, ghost_ids, Chi::mpi.comm);
}

} // namespace chi_physics