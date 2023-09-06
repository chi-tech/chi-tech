#include "MaxMinAvgNodalValuePostProcessor.h"

#include "ChiObjectFactory.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/spatial_discretization.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "mesh/LogicalVolume/LogicalVolume.h"

namespace chi
{

RegisterChiObject(chi, MaxMinAvgNodalValuePostProcessor);

InputParameters MaxMinAvgNodalValuePostProcessor::GetInputParameters()
{
  InputParameters params = PostProcessor::GetInputParameters();
  params += chi_physics::GridBasedFieldFunctionInterface::GetInputParameters();
  params += chi_mesh::LogicalVolumeInterface::GetInputParameters();

  params.SetGeneralDescription(
    "Gets the max/min/avg nodal value of a field function "
    "among nodal values.");
  params.SetDocGroup("doc_PostProcessors");

  params.AddRequiredParameter<std::string>(
    "operation", "The required operation to be performed.");

  using namespace chi_data_types;
  params.ConstrainParameterRange(
    "operation", AllowableRangeList::New({"max", "min", "avg"}));

  return params;
}

MaxMinAvgNodalValuePostProcessor::MaxMinAvgNodalValuePostProcessor(
  const InputParameters& params)
  : PostProcessor(params, PPType::SCALAR),
    chi_physics::GridBasedFieldFunctionInterface(params),
    chi_mesh::LogicalVolumeInterface(params),
    operation_(params.GetParamValue<std::string>("operation"))
{
}

// ##################################################################
void MaxMinAvgNodalValuePostProcessor::Initialize()
{
  const auto* grid_field_function = GetGridBasedFieldFunction();

  ChiLogicalErrorIf(not grid_field_function,
                    "Attempted to access invalid field"
                    "function");

  const auto& grid = grid_field_function->SDM().Grid();

  const auto* logical_volume_ptr_ = GetLogicalVolume();
  if (logical_volume_ptr_ == nullptr)
  {
    cell_local_ids_.reserve(grid.local_cells.size());
    for (const auto& cell : grid.local_cells)
      cell_local_ids_.push_back(cell.local_id_);
  }
  else
  {
    for (const auto& cell : grid.local_cells)
      if (logical_volume_ptr_->Inside(cell.centroid_))
        cell_local_ids_.push_back(cell.local_id_);
  }

  initialized_ = true;
}

// ##################################################################
void MaxMinAvgNodalValuePostProcessor::Execute(const Event& event_context)
{
  if (not initialized_) Initialize();

  const auto* grid_field_function = GetGridBasedFieldFunction();

  ChiLogicalErrorIf(not grid_field_function,
                    "Attempted to access invalid field"
                    "function");

  const auto& ref_ff = *grid_field_function;
  const auto& sdm = ref_ff.SDM();
  const auto& grid = sdm.Grid();

  const auto& uk_man = ref_ff.UnkManager();
  const auto uid = 0;
  const auto cid = 0;

  const auto field_data = ref_ff.GetGhostedFieldVector();
  const size_t num_local_dofs =
    ref_ff.SDM().GetNumLocalDOFs(ref_ff.UnkManager());

  double local_max_value = 0.0;
  double local_min_value = 0.0;
  double local_accumulation = 0.0;
  bool first_local = true;
  for (const uint64_t cell_local_id : cell_local_ids_)
  {
    const auto& cell = grid.local_cells[cell_local_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t imap = sdm.MapDOFLocal(cell, i, uk_man, uid, cid);
      if (imap >= 0 and imap < num_local_dofs)
      {
        const double field_value = field_data[imap];
        if (first_local)
        {
          local_max_value = field_value;
          local_min_value = field_value;
          first_local = false;
        }

        local_max_value = std::max(local_max_value, field_value);
        local_min_value = std::min(local_min_value, field_value);
        local_accumulation += field_value;
      }
    } // for i
  }   // for cell-id

  if (operation_ == "max")
  {
    double globl_max_value;
    MPI_Allreduce(&local_max_value, // sendbuf
                  &globl_max_value, // recvbuf
                  1,
                  MPI_DOUBLE,     // count + datatype
                  MPI_MAX,        // op
                  Chi::mpi.comm); // communicator

    value_ = ParameterBlock("", globl_max_value);
  }
  else if (operation_ == "min")
  {
    double globl_min_value;
    MPI_Allreduce(&local_min_value, // sendbuf
                  &globl_min_value, // recvbuf
                  1,
                  MPI_DOUBLE,     // count + datatype
                  MPI_MIN,        // op
                  Chi::mpi.comm); // communicator

    value_ = ParameterBlock("", globl_min_value);
  }
  else if (operation_ == "avg")
  {
    double globl_accumulation;
    MPI_Allreduce(&local_accumulation, // sendbuf
                  &globl_accumulation, // recvbuf
                  1,
                  MPI_DOUBLE,     // count + datatype
                  MPI_SUM,        // op
                  Chi::mpi.comm); // communicator

    const size_t num_globl_dofs =
      ref_ff.SDM().GetNumGlobalDOFs(ref_ff.UnkManager());
    value_ = ParameterBlock("", globl_accumulation / double(num_globl_dofs));
  }
  else
    ChiLogicalError("Unsupported operation type \"" + operation_ + "\".");
}

} // namespace chi