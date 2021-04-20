#include "chi_unpartitioned_mesh.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

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
    << "Done reading VTU file: \""
    << options.file_name << "\".";

  //======================================== Remove duplicate vertices
  vtkSmartPointer<vtkCleanUnstructuredGrid> cleaner =
    vtkSmartPointer<vtkCleanUnstructuredGrid>::New();
  cleaner->SetInputData(reader->GetOutput());
  cleaner->Update();
  auto ugrid = cleaner->GetOutput();
//  auto ugrid = reader->GetOutput();
  uint64_t total_cell_count  = ugrid->GetNumberOfCells();
  uint64_t total_point_count = ugrid->GetNumberOfPoints();

  chi_log.Log(LOG_0)
    << "Clean grid num cells and points: "
    << total_cell_count << " "
    << total_point_count;

  //======================================== Push cells
  size_t num_polyhedrons  = 0;
  size_t num_hexahedrons  = 0;
  size_t num_tetrahedrons = 0;
  size_t num_polygons     = 0;
  size_t num_quads        = 0;
  size_t num_triangles    = 0;
  for (int c=0; c<total_cell_count; ++c)
  {
    auto vtk_cell = ugrid->GetCell(c);
    auto vtk_celltype = vtk_cell->GetCellType();
    if (vtk_celltype == VTK_POLYHEDRON)
    {
      raw_cells.push_back(CreateCellFromVTKPolyhedron(vtk_cell));
      ++num_polyhedrons;
    }
    else if (vtk_celltype == VTK_HEXAHEDRON)
    {
      raw_cells.push_back(CreateCellFromVTKHexahedron(vtk_cell));
      ++num_hexahedrons;
    }
    else if (vtk_celltype == VTK_TETRA)
    {
      raw_cells.push_back(CreateCellFromVTKTetrahedron(vtk_cell));
      ++num_tetrahedrons;
    }
    else if (vtk_celltype == VTK_POLYGON)
    {
      raw_cells.push_back(CreateCellFromVTKPolygon(vtk_cell));
      ++num_polygons;
    }
    else if (vtk_celltype == VTK_QUAD)
    {
      raw_cells.push_back(CreateCellFromVTKQuad(vtk_cell));
      ++num_quads;
    }
    else if (vtk_celltype == VTK_TRIANGLE)
    {
      raw_cells.push_back(CreateCellFromVTKTriangle(vtk_cell));
      ++num_triangles;
    }
    else
      throw std::invalid_argument(std::string(__FUNCTION__) +
                                  ": Unsupported cell type.");
  }//for c
  chi_log.Log() << "Number cells read: " << total_cell_count << "\n"
    << "polyhedrons  : " << num_polyhedrons  << "\n"
    << "hexahedrons  : " << num_hexahedrons  << "\n"
    << "tetrahedrons : " << num_tetrahedrons << "\n"
    << "polygons     : " << num_polygons     << "\n"
    << "quads        : " << num_quads        << "\n"
    << "triangles    : " << num_triangles;

  //======================================== Push points
  for (int p=0; p<total_point_count; ++p)
  {
    auto point = ugrid->GetPoint(p);
    vertices.push_back(new chi_mesh::Vertex(point[0],point[1],point[2]));

    if (point[0] < bound_box.xmin) bound_box.xmin = point[0];
    if (point[0] > bound_box.xmax) bound_box.xmax = point[0];
    if (point[1] < bound_box.ymin) bound_box.ymin = point[1];
    if (point[1] > bound_box.ymax) bound_box.ymax = point[1];
    if (point[2] < bound_box.zmin) bound_box.zmin = point[2];
    if (point[2] > bound_box.zmax) bound_box.zmax = point[2];

  }

  //======================================== Always do this
  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();
}