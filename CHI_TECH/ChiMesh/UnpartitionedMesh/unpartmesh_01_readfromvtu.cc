#include "chi_unpartitioned_mesh.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog chi_log;
extern ChiMPI chi_mpi;

#include <fstream>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include "vtkCleanUnstructuredGrid.h"
#include <vtkCellType.h>

#include <vtkDataObject.h>
//###################################################################
/**Reads a VTK unstructured mesh.*/
void chi_mesh::UnpartitionedMesh::
  ReadFromVTU(const chi_mesh::UnpartitionedMesh::Options &options)
{
  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);

  if (!file.is_open())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Failed to open file: "<< options.file_name <<" in call "
      << "to ReadFromVTU \n";
    exit(EXIT_FAILURE);
  }
  file.close();

  //======================================== Read the file
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(options.file_name.c_str());
  chi_log.Log(LOG_0)
    << "Reading VTU file     : \""
    << options.file_name << "\".";

  reader->Update();

  chi_log.Log(LOG_0)
    << "Done reading Ensight-Gold file: \""
    << options.file_name << "\".";

  //======================================== Remove duplicate vertices
  vtkSmartPointer<vtkCleanUnstructuredGrid> cleaner =
    vtkSmartPointer<vtkCleanUnstructuredGrid>::New();
  cleaner->SetInputData(reader->GetOutput());
  cleaner->Update();
  auto ugrid = cleaner->GetOutput();
  int total_cell_count  = ugrid->GetNumberOfCells();
  int total_point_count = ugrid->GetNumberOfPoints();

  chi_log.Log(LOG_0)
    << "Clean grid num cells and points: "
    << total_cell_count << " "
    << total_point_count;

  //======================================== Push cells
  for (int c=0; c<total_cell_count; ++c)
  {
    auto vtk_cell = ugrid->GetCell(c);
    auto vtk_celltype = vtk_cell->GetCellType();
    if (vtk_celltype == VTK_POLYHEDRON)
      raw_cells.push_back(CreateCellFromVTKPolyhedron(vtk_cell));
    else if (vtk_celltype == VTK_HEXAHEDRON)
      raw_cells.push_back(CreateCellFromVTKHexahedron(vtk_cell));
    else if (vtk_celltype == VTK_TETRA)
      raw_cells.push_back(CreateCellFromVTKTetrahedron(vtk_cell));
  }//for c

  //======================================== Push points
  for (int p=0; p<total_point_count; ++p)
  {
    auto point = ugrid->GetPoint(p);
    vertices.push_back(new chi_mesh::Vertex(point[0],point[1],point[2]));
  }


}