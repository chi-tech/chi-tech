#include "chi_unpartitioned_mesh.h"

#include <vtkUnstructuredGrid.h>
#include <vtkAppendFilter.h>

#include <vtkCellData.h>

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Copies the vtk data structures to the current object's internal
 * data.*/
void chi_mesh::UnpartitionedMesh::
  CopyUGridCellsAndPoints(vtkUnstructuredGrid& ugrid,
                          const double scale)
{
  const std::string fname =
    "chi_mesh::UnpartitionedMesh::CopyUGridCellsAndPoints";

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
    auto vtk_celldim = vtk_cell->GetCellDimension();

    if (vtk_celldim == 3)
      raw_cells.push_back(CreateCellFromVTKPolyhedron(vtk_cell));
    else if (vtk_celldim == 2)
      raw_cells.push_back(CreateCellFromVTKPolygon(vtk_cell));
    else if (vtk_celldim == 1)
      raw_cells.push_back(CreateCellFromVTKLine(vtk_cell));
    else if (vtk_celldim == 0)
      raw_cells.push_back(CreateCellFromVTKVertex(vtk_cell));
    else
      throw std::logic_error(fname + ": Unsupported cell dimension.");
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


//###################################################################
/**Set material-ids based on block-wise material ids.*/
void chi_mesh::UnpartitionedMesh::
  SetMaterialIDsFromBlocks(const std::vector<uint64_t> &block_mat_ids)
{
  size_t total_cell_count = raw_cells.size();
  for (size_t c=0; c<total_cell_count; ++c)
  {
    auto cell = raw_cells[c];

    int mat_id=-1;
    uint64_t prev_block_lim = 0;
    for (auto block_lim : block_mat_ids)
    {
      ++mat_id;
      if (c<block_lim and c>=prev_block_lim) break;

      prev_block_lim=block_lim;
    }

    cell->material_id = mat_id;
  }
}


//###################################################################
/**Set material-ids from list.*/
void chi_mesh::UnpartitionedMesh::
  SetMaterialIDsFromList(const std::vector<int> &material_ids)
{
  const size_t total_cell_count = raw_cells.size();

  for (size_t c = 0; c < total_cell_count; ++c)
    raw_cells[c]->material_id = material_ids[c];
}


