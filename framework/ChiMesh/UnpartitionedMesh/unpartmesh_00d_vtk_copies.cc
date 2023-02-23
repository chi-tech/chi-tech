#include "chi_unpartitioned_mesh.h"

#include <vtkUnstructuredGrid.h>
#include <vtkAppendFilter.h>
#include <vtkCleanUnstructuredGrid.h>

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/***/
int chi_mesh::UnpartitionedMesh::
  FindHighestDimension(std::vector<vtkUGridPtr> &ugrid_blocks)
{

  int max_dim = 0;
  for (auto& ugrid : ugrid_blocks)
    if (ugrid->GetNumberOfCells() != 0)
      max_dim = std::max(max_dim, ugrid->GetCell(0)->GetCellDimension());

  return max_dim;
}


//###################################################################
/***/
std::vector<uint64_t> chi_mesh::UnpartitionedMesh::
  BuildBlockCellExtents(std::vector<vtkUGridPtr> &ugrid_blocks,
                        const int desired_dimension)
{
  std::vector<uint64_t> block_mat_ids;
  size_t total_cells = 0;

  for (auto& ugrid : ugrid_blocks)
  {
    uint64_t num_cells  = ugrid->GetNumberOfCells();

    if (num_cells == 0) continue;

    if (ugrid->GetCell(0)->GetCellDimension() == desired_dimension)
    {
      total_cells += num_cells;
      block_mat_ids.push_back(total_cells);
    }
  }
  return block_mat_ids;
}


//###################################################################
/***/
chi_mesh::UnpartitionedMesh::vtkUGridPtr chi_mesh::UnpartitionedMesh::
  ConsolidateAndCleanBlocks(std::vector<vtkUGridPtr> &ugrid_blocks,
                            const int desired_dimension)
{
  auto append = vtkSmartPointer<vtkAppendFilter>::New();
  for (auto& ugrid : ugrid_blocks)
  {
    if (ugrid->GetNumberOfCells() == 0) continue;

    if (ugrid->GetCell(0)->GetCellDimension() == desired_dimension)
    {
      append->AddInputData(ugrid);
      append->Update();
    }
  }
  chi::log.Log0Verbose1() << "Updating appended filter.";

  chi::log.Log0Verbose1() << "Getting dirty grid.";
  auto dirty_ugrid = vtkSmartPointer<vtkUnstructuredGrid>(
    vtkUnstructuredGrid::SafeDownCast(append->GetOutput()));

  chi::log.Log0Verbose1()
    << "Dirty grid num cells and points: "
    << dirty_ugrid->GetNumberOfCells() << " "
    << dirty_ugrid->GetNumberOfPoints();

  //======================================== Remove duplicate vertices
  auto cleaner = vtkSmartPointer<vtkCleanUnstructuredGrid>::New();
  cleaner->SetInputData(dirty_ugrid);
  cleaner->Update();
  auto ugrid = vtkUGridPtr(cleaner->GetOutput());

  return ugrid;
}


//###################################################################
/***/
void chi_mesh::UnpartitionedMesh::
  CopyUGridCellsAndPoints(vtkUnstructuredGrid& ugrid,
                          const std::vector<uint64_t>& block_mat_ids,
                          const double scale)
{
  uint64_t total_cell_count  = ugrid.GetNumberOfCells();
  uint64_t total_point_count = ugrid.GetNumberOfPoints();

  chi::log.Log()
    << "Clean grid num cells and points: "
    << total_cell_count << " "
    << total_point_count;

  //======================================== Push cells
  for (size_t c=0; c<total_cell_count; ++c)
  {
    auto vtk_cell = ugrid.GetCell(static_cast<vtkIdType>(c));
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
    for (auto block_lim : block_mat_ids)
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
    auto point = ugrid.GetPoint(static_cast<vtkIdType>(p));

    point[0] = point[0]*scale;
    point[1] = point[1]*scale;
    point[2] = point[2]*scale;

    vertices.emplace_back(point[0],point[1],point[2]);

    if (point[0] < bound_box.xmin) bound_box.xmin = point[0];
    if (point[0] > bound_box.xmax) bound_box.xmax = point[0];
    if (point[1] < bound_box.ymin) bound_box.ymin = point[1];
    if (point[1] > bound_box.ymax) bound_box.ymax = point[1];
    if (point[2] < bound_box.zmin) bound_box.zmin = point[2];
    if (point[2] > bound_box.zmax) bound_box.zmax = point[2];
  }
}