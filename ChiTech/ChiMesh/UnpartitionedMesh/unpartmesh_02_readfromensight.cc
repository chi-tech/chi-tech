#include "chi_unpartitioned_mesh.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

#include <sstream>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkPolyhedron.h>
#include <vtkEnSightGoldBinaryReader.h>

#include <vtkMultiBlockDataSet.h>

#include <vtkAppendFilter.h>

#include <vtkCleanUnstructuredGrid.h>

//###################################################################
/**Reads an Ensight-Gold unstructured mesh.*/
void chi_mesh::UnpartitionedMesh::
  ReadFromEnsightGold(const chi_mesh::UnpartitionedMesh::Options &options)
{
  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);

  if (!file.is_open())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Failed to open file: "<< options.file_name <<" in call "
      << "to ReadFromEnsightGold \n";
    exit(EXIT_FAILURE);
  }
  file.close();

  //======================================== Read the file
  mesh_options = options;
  auto reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
  reader->SetCaseFileName(options.file_name.c_str());
  chi_log.Log(LOG_0)
    << "Reading Ensight-Gold file     : \""
    << options.file_name << "\".";

  reader->Update();

  chi_log.Log(LOG_0)
    << "Done reading Ensight-Gold file: \""
    << options.file_name << "\".";

  //======================================== Separate the blocks
  auto multiblock = reader->GetOutput();
  size_t num_blocks = multiblock->GetNumberOfBlocks();
  chi_log.Log(LOG_0) << "Number of blocks in file: " << num_blocks;

  std::vector<vtkSmartPointer<vtkUnstructuredGrid>> grid_blocks;
  grid_blocks.reserve(num_blocks);
  for (size_t b=0; b<num_blocks; ++b)
    grid_blocks.emplace_back(
        vtkUnstructuredGrid::SafeDownCast(multiblock->GetBlock(b))
      );

  //========================================= Determine mesh type
  // This step goes through each grid block
  // and tries to find 3D cells. If it does the
  // overall mesh is classified as 3D.
  bool mesh_is_3D = false;
  for (auto& ugrid : grid_blocks)
  {
    if (ugrid->GetNumberOfCells() == 0) continue;

    auto cell = ugrid->GetCell(0); //Get first cell
    auto ctype = cell->GetCellType();

    if ((ctype == VTK_POLYHEDRON) or
        (ctype == VTK_HEXAHEDRON) or
        (ctype == VTK_TETRA) )
    {
      mesh_is_3D = true;
      break;
    }
  }//for grid block

  //========================================= Process each block
  auto append = vtkSmartPointer<vtkAppendFilter>::New();
  std::vector<uint64_t> block_mat_id;
  size_t total_cells = 0;
  int ug=-1;
  for (auto& ugrid : grid_blocks)
  {
    ++ug;
    uint64_t num_cells  = ugrid->GetNumberOfCells();
    uint64_t num_points = ugrid->GetNumberOfPoints();

    std::stringstream outstr;
    outstr
      << "Block " << ug << " has "
      << num_cells << " cells and "
      << num_points << " points";

    if (num_cells == 0) continue;

    auto cell = ugrid->GetCell(0);
    auto ctype = cell->GetCellType();

    if (mesh_is_3D)
    {
      if (not ((ctype == VTK_POLYHEDRON) or
               (ctype == VTK_HEXAHEDRON) or
               (ctype == VTK_TETRA)) )
      { outstr << " (Boundary Mesh)"; }
      else
      {
        append->AddInputData(ugrid);
        append->Update();
        total_cells += num_cells;
        block_mat_id.push_back(total_cells);
      }
    }//3D mesh
    else
    {
      if (not ((ctype == VTK_POLYGON) or
               (ctype == VTK_QUAD) or
               (ctype == VTK_TRIANGLE)) )
      { outstr << " (Boundary Mesh)"; }
      else
      {
        append->AddInputData(ugrid);
        append->Update();
        total_cells += num_cells;
        block_mat_id.push_back(total_cells);
      }
    }//2D mesh

   chi_log.Log(LOG_0VERBOSE_1) << outstr.str();
  }
  chi_log.Log(LOG_0VERBOSE_1) << "Updating appended filter.";
//  append->Update();
  chi_log.Log(LOG_0VERBOSE_1) << "Getting dirty grid.";
  auto dirty_ugrid = vtkSmartPointer<vtkUnstructuredGrid>(
    vtkUnstructuredGrid::SafeDownCast(append->GetOutput()));

  chi_log.Log(LOG_0VERBOSE_1)
    << "Dirty grid num cells and points: "
    << dirty_ugrid->GetNumberOfCells() << " "
    << dirty_ugrid->GetNumberOfPoints();

  //======================================== Remove duplicate vertices
  auto cleaner = vtkSmartPointer<vtkCleanUnstructuredGrid>::New();
  cleaner->SetInputData(dirty_ugrid);
  cleaner->Update();
  auto ugrid = cleaner->GetOutput();
  uint64_t total_cell_count  = ugrid->GetNumberOfCells();
  uint64_t total_point_count = ugrid->GetNumberOfPoints();

  chi_log.Log(LOG_0)
    << "Clean grid num cells and points: "
    << total_cell_count << " "
    << total_point_count;

  //======================================== Push cells
  for (size_t c=0; c<total_cell_count; ++c)
  {
    auto vtk_cell = ugrid->GetCell(static_cast<vtkIdType>(c));
    auto vtk_celltype = vtk_cell->GetCellType();
    if (vtk_celltype == VTK_POLYHEDRON)
      raw_cells.push_back(CreateCellFromVTKPolyhedron(vtk_cell));
    else if (vtk_celltype == VTK_HEXAHEDRON)
      raw_cells.push_back(CreateCellFromVTKHexahedron(vtk_cell));
    else if (vtk_celltype == VTK_TETRA)
      raw_cells.push_back(CreateCellFromVTKTetrahedron(vtk_cell));
    else if (vtk_celltype == VTK_POLYGON)
      raw_cells.push_back(CreateCellFromVTKPolygon(vtk_cell));
    else if (vtk_celltype == VTK_QUAD)
      raw_cells.push_back(CreateCellFromVTKQuad(vtk_cell));
    else if (vtk_celltype == VTK_TRIANGLE)
      raw_cells.push_back(CreateCellFromVTKTriangle(vtk_cell));

    int mat_id=-1;
    uint64_t prev_block_lim = 0;
    for (auto block_lim : block_mat_id)
    {
      ++mat_id;
      if (c<block_lim and c>=prev_block_lim) break;

      prev_block_lim=block_lim;
    }

    raw_cells.back()->material_id = mat_id;
  }//for c

  //======================================== Push points
  for (size_t p=0; p<total_point_count; ++p)
  {
    auto point = ugrid->GetPoint(static_cast<vtkIdType>(p));

    point[0] = point[0]*options.scale;
    point[1] = point[1]*options.scale;
    point[2] = point[2]*options.scale;


    vertices.emplace_back(point[0],point[1],point[2]);

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