#include "chi_unpartitioned_mesh.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>


#define ErrorReadingFile(fname) \
std::runtime_error("Failed to open file: " + options.file_name + \
" in call to " + #fname + ".")

//###################################################################
/**Reads a VTK unstructured mesh. This reader will use the following
 * options:
 * - `file_name`, of course.
 * - `material_id_fieldname`, cell data for material_id.*/
void chi_mesh::UnpartitionedMesh::
  ReadFromVTU(const chi_mesh::UnpartitionedMesh::Options &options)
{
  chi::log.Log() << "Reading VTU file: " << options.file_name << ".";

  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open()) throw ErrorReadingFile(ReadFromVTU);
  file.close();

  //======================================== Read the file
  mesh_options = options;
  auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read vtu file");
  reader->Update();
  reader->UpdateInformation();

  //======================================== Separate the blocks
  // For vtu files this is very simple. The
  // output of the reader is an UnstructuredGrid.
  auto ugrid_main = vtkUGridPtr(reader->GetOutput());
  std::vector<vtkSmartPointer<vtkUnstructuredGrid>> grid_blocks = {ugrid_main};

  //======================================== Process blocks
  const int max_dimension = FindHighestDimension(grid_blocks);
  auto ugrid = ConsolidateAndCleanBlocks(grid_blocks, max_dimension);

  //======================================== Copy Data
  CopyUGridCellsAndPoints(*ugrid, options.scale);

  //======================================== Set material ids
  const auto material_ids = BuildCellMaterialIDsFromField(
    ugrid, options.material_id_fieldname, options.file_name);
  SetMaterialIDsFromList(material_ids);

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

  chi::log.Log() << "Done reading VTU file: " << options.file_name << ".";
}

