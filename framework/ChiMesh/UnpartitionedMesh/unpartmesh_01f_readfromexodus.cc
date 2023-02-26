#include "chi_unpartitioned_mesh.h"

#include "ChiMesh/MeshContinuum/chi_grid_vtk_utils.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkExodusIIReader.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkInformation.h>

#include "chi_runtime.h"
#include "chi_log.h"

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
  mesh_options_ = options;
  auto reader = vtkSmartPointer<vtkExodusIIReader>::New();
  reader->SetFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read file-type with this routine");

  reader->UpdateInformation();
  //Exodus ships boundary-ids via SideSets. This allows
  //it to be read from the file
  reader->SetAllArrayStatus(reader->SIDE_SET, 1);
  reader->SetAllArrayStatus(reader->SIDE_SET_CONN, 1);
  reader->Update();

  //======================================== Separate the blocks
  // This part was quite difficult. I eventually found how to do
  // this from this post:
  // https://public.kitware.com/pipermail/vtkusers/2010-December/064921.html
  // It indicated that all exodus formats are read with a single
  // top level vtkMultiBlockDataSet (Level 1). Each of the blocks in
  // level 2 is also of type vtkMultiBlockDataSet. The level 2block that has
  // the actual elements is also split into blocks but these, level 3,
  // blocks each contain a structure castable to vtkUnstructuredGrid.
  auto multiblock = reader->GetOutput();

  std::vector<vtkUGridPtrAndName> grid_blocks;
  auto iter_a = multiblock->NewIterator();
  iter_a->GoToFirstItem();
  while (not iter_a->IsDoneWithTraversal())
  {
    auto block_a = iter_a->GetCurrentDataObject();

    if (block_a->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
      grid_blocks.emplace_back(
        vtkUnstructuredGrid::SafeDownCast(block_a),
        iter_a->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME()));

    iter_a->GoToNextItem();
  }

  //======================================== Process blocks
  const int max_dimension = chi_mesh::FindHighestDimension(grid_blocks);
  auto ugrid = chi_mesh::ConsolidateAndCleanBlocks(grid_blocks, max_dimension);

  //======================================== Copy Data
  CopyUGridCellsAndPoints(*ugrid, options.scale);

  //======================================== Set material ids
  const auto block_mat_ids =
    chi_mesh::BuildBlockCellExtents(grid_blocks, max_dimension);
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

  attributes_ = dimension | UNSTRUCTURED;

  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();

  chi::log.Log() << "Done reading Exodus file: "
                 << options.file_name << ".";
}