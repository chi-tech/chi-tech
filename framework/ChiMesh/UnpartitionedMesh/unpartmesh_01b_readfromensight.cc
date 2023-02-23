#include "chi_unpartitioned_mesh.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkEnSightGoldBinaryReader.h>
#include <vtkMultiBlockDataSet.h>

//###################################################################
/**Reads an Ensight-Gold unstructured mesh.*/
void chi_mesh::UnpartitionedMesh::
  ReadFromEnsightGold(const chi_mesh::UnpartitionedMesh::Options &options)
{
  chi::log.Log() << "Reading Ensight-Gold file: " << options.file_name << ".";

  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);

  if (!file.is_open())
  {
    chi::log.LogAllError()
      << "Failed to open file: "<< options.file_name <<" in call "
      << "to ReadFromEnsightGold \n";
    chi::Exit(EXIT_FAILURE);
  }
  file.close();

  //======================================== Read the file
  mesh_options = options;
  auto reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
  reader->SetCaseFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read ensight file");
  reader->Update();
  reader->UpdateInformation();

  //======================================== Separate the blocks
  auto multiblock = reader->GetOutput();
  size_t num_blocks = multiblock->GetNumberOfBlocks();

  std::vector<vtkSmartPointer<vtkUnstructuredGrid>> grid_blocks;
  grid_blocks.reserve(num_blocks);
  for (size_t b=0; b<num_blocks; ++b)
    grid_blocks.emplace_back(
        vtkUnstructuredGrid::SafeDownCast(multiblock->GetBlock(b)));

  //======================================== Process blocks
  const int max_dimension = FindHighestDimension(grid_blocks);
  const auto block_mat_ids = BuildBlockCellExtents(grid_blocks, max_dimension);
  auto ugrid = ConsolidateAndCleanBlocks(grid_blocks, max_dimension);

  //======================================== Copy Data
  CopyUGridCellsAndPoints(*ugrid, block_mat_ids, options.scale);

  //======================================== Always do this
  chi_mesh::MeshAttributes dimension = NONE;
  switch (max_dimension)
  {
    case 1: dimension = DIMENSION_1; break;
    case 2: dimension = DIMENSION_2; break;
    case 3: dimension = DIMENSION_3; break;
    default: break;
  }

  attributes = dimension | UNSTRUCTURED;

  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();

  chi::log.Log() << "Done reading Ensight-Gold file: "
                 << options.file_name << ".";
}