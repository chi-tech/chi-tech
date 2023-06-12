#include "chi_grid_vtk_utils.h"

#include "chi_meshcontinuum.h"

#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedIntArray.h>

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include "chi_runtime.h"

//###################################################################
/**Uploads vertices and cells to an unstructured grid. This routine
 * also uploads cell material ids (sub-domain ids) and partition ids.*/
vtkNew<vtkUnstructuredGrid> chi_mesh::
  PrepareVtkUnstructuredGrid(const chi_mesh::MeshContinuum& grid)
{
  //============================================= Instantiate VTK items
  vtkNew<vtkUnstructuredGrid>         ugrid;
  vtkNew<vtkPoints>                   points;
  vtkNew<vtkIntArray>                 material_array;
  vtkNew<vtkUnsignedIntArray>         partition_id_array;

  points->SetDataType(VTK_DOUBLE);

  //============================================= Set names
  material_array->SetName("Material");
  partition_id_array->SetName("Partition");

  //############################################# Populate cell information
  int64_t node_count=0;
  for (const auto& cell : grid.local_cells)
  {
    chi_mesh::UploadCellGeometryDiscontinuous(grid, cell, node_count, points, ugrid);

    material_array->InsertNextValue(cell.material_id_);
    partition_id_array->InsertNextValue(cell.partition_id_);
  }//for local cells
  ugrid->SetPoints(points);

  ugrid->GetCellData()->AddArray(material_array);
  ugrid->GetCellData()->AddArray(partition_id_array);

  return ugrid;
}

//###################################################################
/**Writes an unstructured grid to files (.pvtu and .vtu).*/
void chi_mesh::WritePVTUFiles(vtkNew<vtkUnstructuredGrid> &ugrid,
                              const std::string& file_base_name)
{
  //============================================= Construct file name
  std::string base_filename     = std::string(file_base_name);
  std::string location_filename = base_filename +
                                  std::string("_") +
                                  std::to_string(chi::mpi.location_id) +
                                  std::string(".vtu");

  //============================================= Write master file
  if (chi::mpi.location_id == 0)
  {
    std::string pvtu_file_name = base_filename + std::string(".pvtu");

    auto pgrid_writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();

    pgrid_writer->EncodeAppendedDataOff();
    pgrid_writer->SetFileName(pvtu_file_name.c_str());
    pgrid_writer->SetNumberOfPieces(chi::mpi.process_count);
    pgrid_writer->SetStartPiece(chi::mpi.location_id);
    pgrid_writer->SetEndPiece(chi::mpi.process_count-1);
    pgrid_writer->SetInputData(ugrid);

    pgrid_writer->Write();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Serial output each piece
  auto grid_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  grid_writer->Write();
}