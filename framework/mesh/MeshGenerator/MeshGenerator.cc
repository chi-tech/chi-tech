#include "MeshGenerator.h"

#include "mesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"
#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/VolumeMesher/chi_volumemesher.h"
#include "graphs/GraphPartitioner.h"
#include "graphs/PETScGraphPartitioner.h"

#include "ChiObjectFactory.h"

#include "chi_runtime.h"
#include "chi_log.h" // TODO: Remove

namespace chi_mesh
{

RegisterChiObject(chi_mesh, MeshGenerator);

chi::InputParameters MeshGenerator::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription("The base class for all mesh generators");
  params.SetDocGroup("MeshGenerator");

  params.AddOptionalParameter(
    "scale", 1.0, "Uniform scale to apply to the mesh after reading.");

  params.AddOptionalParameterArray(
    "inputs",
    std::vector<size_t>{},
    "A list of handles to MeshGenerator objects");

  params.AddOptionalParameter(
    "partitioner",
    0,
    "Handle to a GraphPartitioner object to use for parallel partitioning."
    "This will default to PETScGraphPartitioner with a \"parmetis\" setting");

  return params;
}

MeshGenerator::MeshGenerator(const chi::InputParameters& params)
  : ChiObject(params),
    scale_(params.GetParamValue<double>("scale"))
{
  //============================================= Convert input handles
  auto input_handles = params.GetParamVectorValue<size_t>("inputs");

  for (const size_t input_handle : input_handles)
  {
    auto& mesh_generator = Chi::GetStackItem<MeshGenerator>(
      Chi::object_stack, input_handle, __FUNCTION__);
    inputs_.push_back(&mesh_generator);
  }

  //============================================= Set partitioner
  size_t partitioner_handle;
  if (params.ParametersAtAssignment().Has("partitioner"))
    partitioner_handle = params.GetParamValue<size_t>("partitioner");
  else
  {
    auto& factory = ChiObjectFactory::GetInstance();
    auto valid_params = chi::PETScGraphPartitioner::GetInputParameters();
    partitioner_handle = factory.MakeRegisteredObjectOfType(
      "chi::PETScGraphPartitioner", chi::ParameterBlock());
  }
  partitioner_ = &Chi::GetStackItem<chi::GraphPartitioner>(
    Chi::object_stack, partitioner_handle, __FUNCTION__);
}

// ##################################################################
/**Default behavior here is to return the input umesh unaltered.*/
std::unique_ptr<UnpartitionedMesh> MeshGenerator::GenerateUnparitionedMesh(
  std::unique_ptr<UnpartitionedMesh> input_umesh)
{
  return input_umesh;
}

/**Final execution step. */
void MeshGenerator::Execute()
{
  //======================================== Execute all input generators
  // Note these could be empty
  std::unique_ptr<UnpartitionedMesh> current_umesh = nullptr;
  for (auto mesh_generator_ptr : inputs_)
  {
    auto new_umesh =
      mesh_generator_ptr->GenerateUnparitionedMesh(std::move(current_umesh));
    current_umesh = std::move(new_umesh);
  }

  //======================================== Generate final umesh and convert it
  current_umesh = GenerateUnparitionedMesh(std::move(current_umesh));
  auto grid_ptr = SetupMesh(std::move(current_umesh));

  //======================================== Assign the mesh to a VolumeMesher
  auto new_mesher =
    std::make_shared<chi_mesh::VolumeMesher>(VolumeMesherType::UNPARTITIONED);
  new_mesher->SetContinuum(grid_ptr);

  if (Chi::current_mesh_handler < 0)
    chi_mesh::PushNewHandlerAndGetIndex();

  auto& cur_hndlr = chi_mesh::GetCurrentHandler();
  cur_hndlr.SetVolumeMesher(new_mesher);

  Chi::mpi.Barrier();
}

} // namespace chi_mesh