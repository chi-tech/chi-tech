#include "fieldfunction2.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>

#include <vtkInformation.h>

//###################################################################
/**Export multiple field functions to VTK.*/
void chi_physics::FieldFunction2::
  ExportMultipleToVTK(
    const std::string &file_base_name,
    const std::vector<std::shared_ptr<const FieldFunction2>> &ff_list)
{
  const std::string fname = "chi_physics::FieldFunction2::ExportMultipleToVTK";
  chi::log.Log() << "Exporting field functions to VTK with file base \""
                 << file_base_name << "\"";

  if (ff_list.empty())
    throw std::logic_error(fname + ": Cannot be used with empty field-function"
                                   " list");

  //============================================= Setup master and check slaves
  const auto& master_ff_ptr = ff_list.front();
  const auto& master_ff = *master_ff_ptr;

  for (const auto& ff_ptr : ff_list)
    if (ff_ptr != master_ff_ptr)
      if (ff_ptr->m_sdm->ref_grid != master_ff_ptr->m_sdm->ref_grid)
        throw std::logic_error(fname +
        ": Cannot be used with field functions based on different grids.");

  //============================================= Get grid
  const auto& grid = *master_ff.m_sdm->ref_grid;

  //============================================= Instantiate VTK items
  vtkNew<vtkUnstructuredGrid>         ugrid;
  vtkNew<vtkPoints>                   points;
  vtkNew<vtkIntArray>                 material_array;
  vtkNew<vtkUnsignedIntArray>         partition_id_array;

  //============================================= Set names
  material_array->SetName("Material");
  partition_id_array->SetName("Partition");

  //############################################# Populate cell information
  int64_t node_count=0;
  for (const auto& cell : grid.local_cells)
  {
    UploadCellGeometry(grid, cell, node_count, points, ugrid);

    material_array->InsertNextValue(cell.material_id);
    partition_id_array->InsertNextValue(cell.partition_id);
  }//for local cells
  ugrid->SetPoints(points);

  ugrid->GetCellData()->AddArray(material_array);
  ugrid->GetCellData()->AddArray(partition_id_array);

  //============================================= Upload cell/point data
  auto cell_data = ugrid->GetCellData();
  auto point_data = ugrid->GetPointData();
  for (const auto& ff_ptr : ff_list)
  {
    const auto field_vector = ff_ptr->GetGhostedFieldVector();

    const auto& uk_man = ff_ptr->m_unknown_manager;
    const auto& unknown = ff_ptr->m_unknown;
    const auto& sdm = ff_ptr->m_sdm;

    for (uint c=0; c<unknown.num_components; ++c)
    {
      const std::string component_name = ff_ptr->m_text_name +
                                         unknown.text_name +
                                         unknown.component_text_names[c];
      vtkNew<vtkDoubleArray> point_array;
      vtkNew<vtkDoubleArray> cell_array;

      point_array->SetName(component_name.c_str());
      cell_array->SetName(component_name.c_str());

      //Populate the array here
      for (const auto& cell : grid.local_cells)
      {
        const size_t num_nodes = sdm->GetCellNumNodes(cell);

        if (num_nodes == cell.vertex_ids.size())
        {
          double node_average = 0.0;
          for (int n=0; n<num_nodes; ++n)
          {
            const int64_t nmap = sdm->MapDOFLocal(cell,n,uk_man,0,c);

            const double field_value = field_vector[nmap];

            point_array->InsertNextValue(field_value);
            node_average += field_value;
          }//for node
          node_average /= static_cast<double>(num_nodes);
          cell_array->InsertNextValue(node_average);
        }
        else
        {
          double node_average = 0.0;
          for (int n=0; n<num_nodes; ++n)
          {
            const int64_t nmap = sdm->MapDOFLocal(cell,n,uk_man,0,c);

            const double field_value = field_vector[nmap];
            node_average += field_value;
          }//for node
          node_average /= static_cast<double>(num_nodes);
          cell_array->InsertNextValue(node_average);
          for (int n=0; n<cell.vertex_ids.size(); ++n)
          {
            point_array->InsertNextValue(node_average);
          }//for vertex
        }

      }//for cell

      point_data->AddArray(point_array);
      cell_data->AddArray(cell_array);
    }//for component
  }//for ff_ptr

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

  chi::log.Log() << "Done exporting field functions to VTK.";
}