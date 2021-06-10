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
    for (int c=0; c < num_components; c++)
    {
      component_names[c] = field_name + ff_uk.component_text_names[c];
      field_node_array[c]->SetName(component_names[c].c_str());
      field_cell_avg_array[c]->SetName((component_names[c] +
                                        std::string("-avg")).c_str());
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
    for (int c=0; c<num_components; ++c)
      cell_comps_to_map[c] = std::make_pair(cell.local_id,(num_components==1)?
                                                          ref_component : c);

    CreateFVMappingLocal(cell_comps_to_map, mapping);

    for (int c=0; c < num_components; ++c)
    {
      double dof_value = field_vector_local->operator[](mapping[c]);
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
  for (int c=0; c<num_components; ++c)
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



////###################################################################
///**Handles the PWLD version of a field function export to VTK.
// *
// * */
//void chi_physics::FieldFunction::ExportToVTKFVG(const std::string& base_name,
//                                                const std::string& field_name)
//{
//  chi_mesh::FieldFunctionInterpolation ff_interpol;
//
//  std::vector<std::vector<double>>    d_nodes;
//
//  auto points = vtkSmartPointer<vtkPoints>::New();
//
//  const auto& field_unknown = unknown_manager.unknowns[ref_variable];
//
//  //============================================= Init grid and material name
//  vtkUnstructuredGrid* ugrid;
//  vtkIntArray*      matarray;
//  vtkIntArray*      pararray;
//
//  std::vector<vtkDoubleArray*>   phiavgarray(field_unknown.num_components);
//
//  ugrid    = vtkUnstructuredGrid::New();
//  matarray = vtkIntArray::New();
//  matarray->SetName("Material");
//  pararray = vtkIntArray::New();
//  pararray->SetName("Partition");
//
//  for (int g=0; g < field_unknown.num_components; g++)
//  {
//    phiavgarray[g] = vtkDoubleArray::New();
//    phiavgarray[g]->SetName((field_name +
//                             std::string("_g") +
//                             std::to_string(g) +
//                             std::string("-Avg")).c_str());
//  }
//
//
//  //========================================= Populate nodes
//  size_t vcount=0;
//  for (const auto& vertex : grid->vertices)
//  {
//    std::vector<double> d_node;
//    d_node.push_back(vertex.x);
//    d_node.push_back(vertex.y);
//    d_node.push_back(vertex.z);
//
//    d_nodes.push_back(d_node);
//
//    points->InsertPoint(static_cast<vtkIdType>(vcount),d_node.data());
//    ++vcount;
//  }
//
//  const auto& ff_uk = this->unknown_manager.unknowns[ref_variable];
//
//  //======================================== populate cell mapping
//  size_t num_loc_cells = grid->local_cells.size();
//  std::vector<std::pair<uint64_t,uint>> cell_comps_to_map(num_loc_cells);
//  std::vector<uint64_t> mapping;
//  for (uint64_t lc=0; lc<num_loc_cells; lc++)
//    for (unsigned int g=0; g < ff_uk.num_components; g++)
//      cell_comps_to_map[lc] = std::make_pair(lc,g);
//
//  CreateFVMappingLocal(cell_comps_to_map, mapping);
//
//  //======================================== Populate cell information
//  int counter=-1;
//  for (const auto& cell : grid->local_cells)
//  {
//    int mat_id = cell.material_id;
//
//    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
//    if (cell.Type() == chi_mesh::CellType::SLAB)
//    {
//      auto slab_cell = (chi_mesh::CellSlab*)(&cell);
//
//      std::vector<vtkIdType> cell_info;
//      cell_info.push_back(static_cast<vtkIdType>(slab_cell->vertex_ids[0]));
//      cell_info.push_back(static_cast<vtkIdType>(slab_cell->vertex_ids[1]));
//
//      ugrid->
//        InsertNextCell(VTK_LINE,2,
//                       cell_info.data());
//
//      matarray->InsertNextValue(mat_id);
//      pararray->InsertNextValue(static_cast<int>(cell.partition_id));
//
//      for (int g=0; g < ff_uk.num_components; g++)
//      {
//        ++counter;
//        double phi_value = field_vector_local->operator[](mapping[counter]);
//        phiavgarray[g]->InsertNextValue(phi_value);
//      }//for g
//
//    }
//
//    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
//    if (cell.Type() == chi_mesh::CellType::POLYGON)
//    {
//      auto poly_cell = (chi_mesh::CellPolygon*)(&cell);
//
//      std::vector<vtkIdType> cell_info;
//
//      size_t num_verts = poly_cell->vertex_ids.size();
//      for (int v=0; v<num_verts; v++)
//        cell_info.push_back(static_cast<vtkIdType>(poly_cell->vertex_ids[v]));
//
//      ugrid->
//        InsertNextCell(VTK_POLYGON,
//                       static_cast<vtkIdType>(num_verts),
//                       cell_info.data());
//
//      matarray->InsertNextValue(mat_id);
//      pararray->InsertNextValue(static_cast<int>(cell.partition_id));
//
//      for (int g=0; g < ff_uk.num_components; g++)
//      {
//        ++counter;
//        double phi_value = field_vector_local->operator[](mapping[counter]);
//        phiavgarray[g]->InsertNextValue(phi_value);
//      }//for g
//    }
//
//    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
//    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
//    {
//      auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);
//
//      size_t num_verts = polyh_cell->vertex_ids.size();
//      std::vector<vtkIdType> cell_info(num_verts);
//      for (int v=0; v<num_verts; v++)
//        cell_info[v] = static_cast<vtkIdType>(polyh_cell->vertex_ids[v]);
//
//      vtkNew<vtkIdList> faces;
//
//      size_t num_faces = polyh_cell->faces.size();
//      for (int f=0; f<num_faces; f++)
//      {
//        size_t num_fverts = polyh_cell->faces[f].vertex_ids.size();
//        std::vector<vtkIdType> face(num_fverts);
//        for (int fv=0; fv<num_fverts; fv++)
//          face[fv] = static_cast<vtkIdType>(polyh_cell->faces[f].vertex_ids[fv]);
//
//        faces->InsertNextId(static_cast<vtkIdType>(num_fverts));
//        for (auto vid : face)
//          faces->InsertNextId(vid);
//      }//for f
//
//      ugrid->
//        InsertNextCell(VTK_POLYHEDRON,
//                       static_cast<vtkIdType>(num_verts),
//                       cell_info.data(),
//                       static_cast<vtkIdType>(num_faces),
//                       faces->GetPointer(0));
//
//      matarray->InsertNextValue(mat_id);
//      pararray->InsertNextValue(static_cast<int>(cell.partition_id));
//
//      for (int g=0; g < ff_uk.num_components; g++)
//      {
//        ++counter;
//        double phi_value = field_vector_local->operator[](mapping[counter]);
//        phiavgarray[g]->InsertNextValue(phi_value);
//      }//for g
//    }//polyhedron
//  }//for local cells
//
//  ugrid->SetPoints(points);
//
//  //============================================= Construct file name
//  std::string base_filename     = std::string(base_name);
//  std::string location_filename = base_filename +
//                                  std::string("_") +
//                                  std::to_string(chi_mpi.location_id) +
//                                  std::string(".vtu");
//
//  //============================================= Serial Output each piece
//  vtkXMLUnstructuredGridWriter* grid_writer =
//    vtkXMLUnstructuredGridWriter::New();
//
//  ugrid->GetCellData()->AddArray(matarray);
//  ugrid->GetCellData()->AddArray(pararray);
//  for (int g=0; g < ff_uk.num_components; g++)
//    ugrid->GetCellData()->AddArray(phiavgarray[g]);
//
//  grid_writer->SetInputData(ugrid);
//  grid_writer->SetFileName(location_filename.c_str());
//
//  /*It seems that cluster systems throw an error when the pvtu file
//   * also tries to write to the serial file.*/
////  if (chi_mpi.location_id != 0)
//    grid_writer->Write();
//
//  //============================================= Parallel summary file
//  if (chi_mpi.location_id == 0)
//  {
//      WritePVTU(base_filename, field_name, (int)ff_uk.num_components);
//  }
//}
