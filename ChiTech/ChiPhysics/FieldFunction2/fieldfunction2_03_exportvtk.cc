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

void chi_physics::FieldFunction2::
  ExportToVTK(const std::string &file_base_name) const
{
  chi::log.Log() << "Exporting field function to VTK with file base \""
                 << file_base_name << "\"";

  const auto& grid = *m_sdm->ref_grid;

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
//  const auto field_vector = GetGhostedFieldVector();
  const auto& field_vector = m_field_vector;
  auto cell_data = ugrid->GetCellData();
  auto point_data = ugrid->GetPointData();
  for (uint c=0; c<m_unknown.num_components; ++c)
  {
    const std::string component_name = m_text_name + m_unknown.text_name +
                                       m_unknown.component_text_names[c];
    vtkNew<vtkDoubleArray> point_array;
    vtkNew<vtkDoubleArray> cell_array;

    point_array->SetName(component_name.c_str());
    cell_array->SetName(component_name.c_str());

    //Populate the array here
    for (const auto& cell : grid.local_cells)
    {
      const size_t num_nodes = m_sdm->GetCellNumNodes(cell);

      if (num_nodes == cell.vertex_ids.size())
      {
        double node_average = 0.0;
        for (int n=0; n<num_nodes; ++n)
        {
          const int64_t nmap = m_sdm->MapDOFLocal(cell,n,m_unknown_manager,0,c);

          const double field_value = field_vector[nmap];

          point_array->InsertNextValue(field_value);
          node_average += field_value;
        }//for node
        node_average /= static_cast<double>(num_nodes);
        cell_array->InsertNextValue(node_average);
      }
      else
      {
        std::cout << "Nope\n";
        double node_average = 0.0;
        for (int n=0; n<num_nodes; ++n)
        {
          const int64_t nmap = m_sdm->MapDOFLocal(cell,n,m_unknown_manager,0,c);

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

  chi::log.Log() << "Done exporting field function to VTK.";
}