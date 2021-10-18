#include "fieldfunction.h"

#include "ChiPhysics/chi_physics.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;
extern ChiPhysics&  chi_physics_handler;

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>

#include <vtkInformation.h>



//###################################################################
/**Handles the PWLD version of a field function export to VTK.
 *
 * */
void chi_physics::FieldFunction::ExportToVTKFV(const std::string& base_name,
                                               const std::string& field_name,
                                               bool all_components/*=false*/)
{
  //============================================= Init vtk items
  const auto& ff_uk = this->unknown_manager.unknowns[ref_variable];

  size_t num_components = all_components? ff_uk.num_components : 1;
  vtkNew<vtkUnstructuredGrid>         ugrid;
  vtkNew<vtkPoints>                   points;
  vtkNew<vtkIntArray>                 material_array;
  vtkNew<vtkUnsignedIntArray>         partition_id_array;
  std::vector<vtkNew<vtkDoubleArray>> field_node_array(num_components);
  std::vector<vtkNew<vtkDoubleArray>> field_cell_avg_array(num_components);

  //============================================= Set names
  material_array->SetName("Material");
  partition_id_array->SetName("Partition");
  std::vector<std::string> component_names(num_components);
  if (not all_components)
  {
    component_names.back() = field_name;
    field_node_array.back()->SetName(field_name.c_str());
    field_cell_avg_array.back()->SetName((field_name +
                                          std::string("-avg")).c_str());
  }
  else
  {
    for (size_t c=0; c < num_components; ++c)
    {
      component_names[c] = field_name + ff_uk.component_text_names[c];
      field_node_array[c]->SetName(component_names[c].c_str());
      field_cell_avg_array[c]->SetName((component_names[c] +
                                        std::string("-avg")).c_str());
    }//for c
  }

  //############################################# Populate cell information
  int64_t node_count=0;
  for (const auto& cell : grid->local_cells)
  {
    UploadCellGeometry(cell, node_count, points, ugrid);

    material_array->InsertNextValue(cell.material_id);
    partition_id_array->InsertNextValue(cell.partition_id);

    //Populate Arrays
    std::vector<std::pair<uint64_t,uint>> cell_comps_to_map(num_components);
    std::vector<uint64_t> mapping;
    for (size_t c=0; c<num_components; ++c)
      cell_comps_to_map[c] = std::make_pair(cell.local_id,(num_components==1)?
                                                          ref_component : c);

    CreateFVMappingLocal(cell_comps_to_map, mapping);

    for (size_t c=0; c < num_components; ++c)
    {
      double dof_value = field_vector_local->operator[](mapping[c]);
      for (size_t v=0; v<cell.vertex_ids.size(); ++v)
        field_node_array[c]->InsertNextValue(dof_value);
      field_cell_avg_array[c]->InsertNextValue(dof_value);
    }//for component
  }//for local cells

  ugrid->SetPoints(points);

  //============================================= Construct file name
  std::string base_filename     = std::string(base_name);
  std::string location_filename = base_filename +
                                  std::string("_") +
                                  std::to_string(chi_mpi.location_id) +
                                  std::string(".vtu");

  //============================================= Serial Output each piece
  vtkNew<vtkXMLUnstructuredGridWriter> grid_writer;

  ugrid->GetCellData()->AddArray(material_array);
  ugrid->GetCellData()->AddArray(partition_id_array);
  for (size_t c=0; c<num_components; ++c)
  {
    ugrid->GetPointData()->AddArray(field_node_array[c]);
    ugrid->GetCellData()->AddArray(field_cell_avg_array[c]);
  }

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  grid_writer->Write();

  //============================================= Parallel summary file
  if (chi_mpi.location_id == 0)
      WritePVTU(base_filename, field_name, component_names);
}
