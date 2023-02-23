#include "chi_unpartitioned_mesh.h"

#include <vtkUnstructuredGrid.h>
#include <vtkAppendFilter.h>
#include <vtkCleanUnstructuredGrid.h>

#include <vtkCellData.h>

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Finds the highest dimension across all the grid blocks. This
 * is useful when a vtk-read mesh contains multiple blocks. Some of which
 * are boundary faces.*/
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
/**Consolidates all blocks containing cells with the desired dimension.
 * Thereafter it removes duplicate vertices.*/
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
    auto vtk_celltype = vtk_cell->GetCellType();
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

//    if (vtk_celltype == VTK_POLYHEDRON)
//      raw_cells.push_back(CreateCellFromVTKPolyhedron(vtk_cell));
//    else if (vtk_celltype == VTK_HEXAHEDRON)
//      raw_cells.push_back(CreateCellFromVTKHexahedron(vtk_cell));
//    else if (vtk_celltype == VTK_TETRA)
//      raw_cells.push_back(CreateCellFromVTKTetrahedron(vtk_cell));
//    else if (vtk_celltype == VTK_POLYGON)
//      raw_cells.push_back(CreateCellFromVTKPolygon(vtk_cell));
//    else if (vtk_celltype == VTK_QUAD)
//      raw_cells.push_back(CreateCellFromVTKQuad(vtk_cell));
//    else if (vtk_celltype == VTK_TRIANGLE)
//      raw_cells.push_back(CreateCellFromVTKTriangle(vtk_cell));
//    else if (vtk_celltype == VTK_LINE)
//      raw_cells.push_back(CreateCellFromVTKLine(vtk_cell));
//    else if (vtk_celltype == VTK_VERTEX)
//      raw_cells.push_back(CreateCellFromVTKVertex(vtk_cell));
//    else
//      throw std::logic_error(fname + ": Unsupported cell type.");
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
/**Given several unstructured grid blocks, each denoting a material id,
 * this function sets material ids accordingly.*/
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
/**Retrieves material-ids from a field.*/
std::vector<int> chi_mesh::UnpartitionedMesh::
  BuildCellMaterialIDsFromField(vtkUGridPtr &ugrid,
                                const std::string& field_name,
                                const std::string& file_name) const
{
  const size_t total_cell_count = raw_cells.size();
  std::vector<int> material_ids(total_cell_count, -1);

  //======================================== Determine if reading
  //                                         cell identifiers
  vtkDataArray* cell_id_array_ptr;
  if (field_name.empty())
  {
    chi::log.Log0Warning()
      << "A user-supplied field name from which to recover cell identifiers "
      << "has not been provided. Only the mesh will be read.";
    goto end_error_checks;
  }
  else
  {
    auto cell_data = ugrid->GetCellData();
    const auto vtk_abstract_array_ptr =
      cell_data->GetAbstractArray(field_name.c_str());

    if (!vtk_abstract_array_ptr)
    {
      chi::log.Log0Warning()
        << "The VTU file : \"" << file_name << "\" "
        << "does not contain a vtkCellData field of name : \""
        << field_name << "\". Material-ids will be left unassigned.";
      goto end_error_checks;
    }

    cell_id_array_ptr = vtkArrayDownCast<vtkDataArray>(vtk_abstract_array_ptr);
    if (!cell_id_array_ptr)
    {
      chi::log.Log0Warning()
        << "The VTU file : \"" << file_name << "\" "
        << "with vtkCellData field of name : \""
        << field_name << "\" "
        << "cannot be downcast to vtkDataArray. Material-ids will be left "
           "unassigned.";
      goto end_error_checks;
    }

    const auto cell_id_n_tup = cell_id_array_ptr->GetNumberOfTuples();
    if (cell_id_n_tup != total_cell_count)
    {
      chi::log.Log0Warning()
        << "The VTU file : \"" << file_name << "\" "
        << "with vtkCellData field of name : \""
        << field_name << "\" has n. tuples : "
        << cell_id_n_tup << ", but differs from the value expected : "
        << total_cell_count << ". Material-ids will be left unassigned.";
      goto end_error_checks;
    }

    const auto cell_id_n_val = cell_id_array_ptr->GetNumberOfValues();
    if (cell_id_n_val != total_cell_count)
    {
      chi::log.Log0Warning()
        << "The VTU file : \"" << file_name << "\" "
        << "with vtkCellData field of name : \""
        << field_name << "\" has n. values : "
        << cell_id_n_val << ", but differs from the value expected : "
        << total_cell_count << ". Material-ids will be left unassigned.";
      goto end_error_checks;
    }
  }

  //  apply cell identifier
  for (size_t c = 0; c < total_cell_count; ++c)
  {
    std::vector<double> cell_id_vec(1);
    cell_id_array_ptr->GetTuple(static_cast<vtkIdType>(c), cell_id_vec.data());
    const auto mat_id = (int) cell_id_vec.front();

    material_ids[c] = mat_id;
  }

end_error_checks:
  return material_ids;
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


