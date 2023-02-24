#include "chi_grid_vtk_utils.h"

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
int chi_mesh::FindHighestDimension(std::vector<vtkUGridPtrAndName> &ugrid_blocks)
{

  int max_dim = 0;
  for (auto& ugrid : ugrid_blocks)
    if (ugrid.first->GetNumberOfCells() != 0)
      max_dim = std::max(max_dim, ugrid.first->GetCell(0)->GetCellDimension());

  return max_dim;
}


//###################################################################
/**Consolidates all blocks containing cells with the desired dimension.
 * Thereafter it removes duplicate vertices.*/
chi_mesh::vtkUGridPtr chi_mesh::
  ConsolidateAndCleanBlocks(std::vector<vtkUGridPtrAndName> &ugrid_blocks,
                            const int desired_dimension)
{
  auto append = vtkSmartPointer<vtkAppendFilter>::New();
  for (auto& ugrid : ugrid_blocks)
  {
    if (ugrid.first->GetNumberOfCells() == 0) continue;

    if (ugrid.first->GetCell(0)->GetCellDimension() == desired_dimension)
    {
      append->AddInputData(ugrid.first);
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
/**Provides a map of the different grids that have the
 * requested dimension.*/
std::map<std::string, chi_mesh::vtkUGridPtr> chi_mesh::
  SeparateBlocks(std::vector<vtkUGridPtrAndName> &ugrid_blocks,
                 int desired_dimension)
{
  for (auto& ugrid : ugrid_blocks)
  {
    if (ugrid.first->GetNumberOfCells() == 0) continue;

    if (ugrid.first->GetCell(0)->GetCellDimension() == desired_dimension)
    {
//      auto grid_name = ugrid->
    }
  }

  return {};
}


//###################################################################
/**Given several unstructured grid blocks, each denoting a material id,
 * this function sets material ids accordingly.*/
std::vector<uint64_t> chi_mesh::
  BuildBlockCellExtents(std::vector<vtkUGridPtrAndName> &ugrid_blocks,
                        const int desired_dimension)
{
  std::vector<uint64_t> block_mat_ids;
  size_t total_cells = 0;

  for (auto& ugrid : ugrid_blocks)
  {
    uint64_t num_cells  = ugrid.first->GetNumberOfCells();

    if (num_cells == 0) continue;

    if (ugrid.first->GetCell(0)->GetCellDimension() == desired_dimension)
    {
      total_cells += num_cells;
      block_mat_ids.push_back(total_cells);
    }
  }
  return block_mat_ids;
}


//###################################################################
/**Retrieves material-ids from a field.*/
std::vector<int> chi_mesh::
  BuildCellMaterialIDsFromField(vtkUGridPtr &ugrid,
                                const std::string& field_name,
                                const std::string& file_name)
{
  const size_t total_cell_count = ugrid->GetNumberOfCells();
  std::vector<int> material_ids(total_cell_count, -1);

  //======================================== Determine if reading
  //                                         cell identifiers
  vtkDataArray* cell_id_array_ptr;
  if (field_name.empty())
  {
    chi::log.Log0Warning()
      << "A user-supplied field name from which to recover material identifiers "
      << "has not been found. Material-ids will be left unassigned.";
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