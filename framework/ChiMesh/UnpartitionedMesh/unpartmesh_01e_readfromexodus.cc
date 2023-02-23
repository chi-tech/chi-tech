#include "chi_unpartitioned_mesh.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkExodusIIReader.h>
#include <vtkMultiBlockDataSet.h>

#define ErrorReadingFile(fname) \
std::runtime_error("Failed to open file: " + options.file_name + \
" in call to " + #fname + ".")

//###################################################################
/**Reads an Exodus unstructured mesh.*/
void chi_mesh::UnpartitionedMesh::
  ReadFromExodus(const chi_mesh::UnpartitionedMesh::Options &options)
{
  chi::log.Log() << "Reading Exodus file: " << options.file_name << ".";

  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open()) throw ErrorReadingFile(ReadFromExodus);
  file.close();

  //======================================== Read the file
  mesh_options = options;
  auto reader = vtkSmartPointer<vtkExodusIIReader>::New();
  reader->SetFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read exodus file");
  reader->Update();
  reader->UpdateInformation();

  //======================================== Separate the blocks
  // This part was quite difficult. I eventually found how to do
  // this from this post:
  // https://public.kitware.com/pipermail/vtkusers/2010-December/064921.html
  // It indicated that all exodus formats are read with a single
  // top level vtkMultiBlockDataSet (Level 1). Each of the blocks in
  // level 2 is also of type vtkMultiBlockDataSet. The level 2block that has
  // the actual elements is also split into blocks but these, level 3,
  // blocks each contain a structure castable to vtkUnstructuredGrid.
  auto multiblock1 = vtkSmartPointer<vtkMultiBlockDataSet>(reader->GetOutput());
  std::vector<vtkSmartPointer<vtkUnstructuredGrid>> grid_blocks;
  for (int a=0; a < multiblock1->GetNumberOfBlocks(); ++a)
  {
    auto multiblock2 = vtkMultiBlockDataSet::SafeDownCast(
      multiblock1->GetBlock(a));

    if (multiblock2->GetNumberOfCells() == 0) continue;

    for (int b=0; b < multiblock2->GetNumberOfBlocks(); ++b)
      grid_blocks.emplace_back(
        vtkUnstructuredGrid::SafeDownCast(multiblock2->GetBlock(b)));
  }

  //======================================== Process blocks
  const int max_dimension = FindHighestDimension(grid_blocks);
  auto ugrid = ConsolidateAndCleanBlocks(grid_blocks, max_dimension);

  //======================================== Copy Data
  CopyUGridCellsAndPoints(*ugrid, options.scale);

  //======================================== Set material ids
  const auto block_mat_ids =
    BuildBlockCellExtents(grid_blocks, max_dimension);
  SetMaterialIDsFromBlocks(block_mat_ids);

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

  chi::log.Log() << "Done reading Exodus file: "
                 << options.file_name << ".";
}