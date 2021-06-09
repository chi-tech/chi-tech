#include "fieldfunction.h"

#include "ChiMesh/Cell/cell.h"
#include "ChiPhysics/chi_physics.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog&     chi_log;
extern ChiMPI&     chi_mpi;
extern ChiPhysics& chi_physics_handler;

#include <vtkCellType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>

#include <vtkInformation.h>

//###################################################################
/**Handles the PWLD version of a field function export to VTK.
 *
 * */
void chi_physics::FieldFunction::ExportToVTKPWLD(const std::string& base_name,
                                                 const std::string& field_name)
{
  if (spatial_discretization->type !=
      chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Field function spatial discretization"
                                " is not of type "
                                " PIECEWISE_LINEAR_DISCONTINUOUS.");


  //============================================= Init vtk items
  vtkNew<vtkUnstructuredGrid> ugrid;
  vtkNew<vtkPoints>           points;
  vtkNew<vtkIntArray>         material_array;
  vtkNew<vtkUnsignedIntArray> partition_id_array;
  vtkNew<vtkDoubleArray>      field_node_array;
  vtkNew<vtkDoubleArray>      field_cell_avg_array;

  //============================================= Set names
  material_array->SetName("Material");
  partition_id_array->SetName("Partition");
  field_node_array->SetName(field_name.c_str());
  field_cell_avg_array->SetName((field_name + std::string("-Avg")).c_str());


  //############################################# Lambda to convert and add
  //                                              vertices to vtk
  auto AddVerticesToPoints = [this,&points](const size_t num_verts,
                                            const std::vector<uint64_t>& vertex_ids,
                                            int64_t& node_counter)
  {
    std::vector<vtkIdType> cell_info(num_verts);
    for (int v=0; v<num_verts; v++)
    {
      uint64_t vgi = vertex_ids[v];
      std::vector<double> d_node(3);
      d_node[0] = grid->vertices[vgi].x;
      d_node[1] = grid->vertices[vgi].y;
      d_node[2] = grid->vertices[vgi].z;


      points->InsertPoint(node_counter,d_node.data());
      cell_info[v] = node_counter++;
    }

    return cell_info;
  };


  //############################################# Lambda to populate field_node_array
  //                                              and field_cell_avg_array
  auto PopulateArrays = [this](const size_t      num_verts,
                               const uint64_t    cell_local_id,
                               vtkDoubleArray*   phiarray,
                               vtkDoubleArray*   phiavgarray)
  {
    std::vector<uint64_t> mapping;
    std::vector<std::tuple<uint64_t,uint,uint>> cell_node_component_tuples;

    for (int v=0; v<num_verts; v++)
      cell_node_component_tuples.emplace_back(cell_local_id,v,0);

    CreatePWLDMappingLocal(cell_node_component_tuples, mapping);

    double cell_avg_value = 0.0;
    for (int v=0; v<num_verts; v++)
    {
      double dof_value = field_vector_local->operator[](mapping[v]);
      cell_avg_value+= dof_value;
      phiarray->InsertNextValue(dof_value);
    }
    phiavgarray->InsertNextValue(cell_avg_value/static_cast<double>(num_verts));
  };


  //############################################# Populate cell information
  int64_t node_count=0;
  for (const auto& cell : grid->local_cells)
  {
    size_t num_verts = cell.vertex_ids.size();
    auto cell_vids = AddVerticesToPoints(num_verts, cell.vertex_ids, node_count);

    material_array->InsertNextValue(cell.material_id);
    partition_id_array->InsertNextValue(cell.partition_id);

    PopulateArrays(num_verts, cell.local_id, field_node_array, field_cell_avg_array);

    //=========================================== Upload cell geometry
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      ugrid->InsertNextCell(VTK_LINE,
                            static_cast<vtkIdType>(num_verts),
                            cell_vids.data());
    }
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      ugrid-> InsertNextCell(VTK_POLYGON,
                             static_cast<vtkIdType>(num_verts),
                             cell_vids.data());
    }
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      // Build polyhedron faces
      std::vector<vtkIdType> faces_vids;

      size_t num_faces = cell.faces.size();
      for (auto& face : cell.faces)
      {
        size_t num_fverts = face.vertex_ids.size();
        std::vector<vtkIdType> face_info(num_fverts);
        for (int fv=0; fv<num_fverts; fv++)
        {
          int v = 0;
          for (int cv=0; cv<num_verts; ++cv)
            if (cell.vertex_ids[cv] == face.vertex_ids[fv])
            { v = cv; break; }

          face_info[fv] = cell_vids[v];
        }

        faces_vids.push_back(static_cast<vtkIdType>(num_fverts));
        for (auto vid : face_info)
          faces_vids.push_back(vid);
      }//for f

      ugrid-> InsertNextCell(VTK_POLYHEDRON,
                             static_cast<vtkIdType>(num_verts),
                             cell_vids.data(),
                             static_cast<vtkIdType>(num_faces),
                             faces_vids.data());
    }//polyhedron
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
  ugrid->GetPointData()->AddArray(field_node_array);
  ugrid->GetCellData()->AddArray(field_cell_avg_array);

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  grid_writer->Write();

  //============================================= Parallel summary file
  if (chi_mpi.location_id == 0)
      WritePVTU(base_filename, field_name);
}



//###################################################################
/**Handles the PWLD version of a field function export to VTK with all groups.
 *
 * */
void chi_physics::FieldFunction::ExportToVTKPWLDG(const std::string& base_name,
                                                  const std::string& field_name)
{
  if (spatial_discretization->type !=
      chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Field function spatial discretization"
                                " is not of type "
                                " PIECEWISE_LINEAR_DISCONTINUOUS.");

  //============================================= Init vtk items
  const auto& ff_uk = this->unknown_manager.unknowns[ref_variable];

  vtkNew<vtkUnstructuredGrid>         ugrid;
  vtkNew<vtkPoints>                   points;
  vtkNew<vtkIntArray>                 material_array;
  vtkNew<vtkUnsignedIntArray>         partition_id_array;
  std::vector<vtkNew<vtkDoubleArray>> field_node_array(ff_uk.num_components);
  std::vector<vtkNew<vtkDoubleArray>> field_cell_avg_array(ff_uk.num_components);

  //============================================= Set names
  material_array->SetName("Material");
  partition_id_array->SetName("Partition");

  for (int g=0; g < ff_uk.num_components; g++)
  {
    char group_text[100];
    sprintf(group_text,"%03d",g);

    field_node_array[g]   ->SetName((field_name +
                                     std::string("_g") +
                                     std::string(group_text)).c_str());
    field_cell_avg_array[g]->SetName((field_name +
                                      std::string("_g") +
                                      std::string(group_text) +
                                      std::string("_avg")).c_str());
  }


  //############################################# Lambda to convert and add
  //                                              vertices to vtk
  auto AddVerticesToPoints = [this,&points](const size_t num_verts,
                                            const std::vector<uint64_t>& vertex_ids,
                                            int64_t& node_counter)
  {
    std::vector<vtkIdType> cell_info(num_verts);
    for (int v=0; v<num_verts; v++)
    {
      uint64_t vgi = vertex_ids[v];
      std::vector<double> d_node(3);
      d_node[0] = grid->vertices[vgi].x;
      d_node[1] = grid->vertices[vgi].y;
      d_node[2] = grid->vertices[vgi].z;


      points->InsertPoint(node_counter,d_node.data());
      cell_info[v] = node_counter++;
    }

    return cell_info;
  };


  //############################################# Lambda to populate field_node_array
  //                                              and field_cell_avg_array
  auto PopulateArrays = [this,&ff_uk](
    const size_t      num_verts,
    const uint64_t    cell_local_id,
    std::vector<vtkNew<vtkDoubleArray>>& field_node_array,
    std::vector<vtkNew<vtkDoubleArray>>& field_cell_avg_array)
  {
    std::vector<uint64_t> mapping;
    std::vector<std::tuple<uint64_t,uint,uint>> cell_node_component_tuples;

    for (int g=0; g < ff_uk.num_components; g++)
      for (int v=0; v<num_verts; v++)
        cell_node_component_tuples.emplace_back(cell_local_id,v,g);

    CreatePWLDMappingLocal(cell_node_component_tuples, mapping);

    int counter=-1;
    for (int g=0; g < ff_uk.num_components; g++)
    {
      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        ++counter;
        double dof_value = field_vector_local->operator[](mapping[counter]);
        cell_avg_value+= dof_value;
        field_node_array[g]->InsertNextValue(dof_value);
      }
      field_cell_avg_array[g]->InsertNextValue(cell_avg_value/static_cast<double>(num_verts));
    }//for g
  };


  //======================================== Populate cell information
  int64_t node_count=0;
  for (const auto& cell : grid->local_cells)
  {
    size_t num_verts = cell.vertex_ids.size();
    auto cell_vids = AddVerticesToPoints(num_verts, cell.vertex_ids, node_count);

    material_array->InsertNextValue(cell.material_id);
    partition_id_array->InsertNextValue(cell.partition_id);

    PopulateArrays(num_verts, cell.local_id, field_node_array, field_cell_avg_array);

    //=========================================== Upload cell geometry
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      ugrid->InsertNextCell(VTK_LINE,
                            static_cast<vtkIdType>(num_verts),
                            cell_vids.data());
    }
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      ugrid->InsertNextCell(VTK_POLYGON,
                       static_cast<vtkIdType>(num_verts),
                       cell_vids.data());
    }
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      // Build polyhedron faces
      std::vector<vtkIdType> faces_vids;

      size_t num_faces = cell.faces.size();
      for (auto& face : cell.faces)
      {
        size_t num_fverts = face.vertex_ids.size();
        std::vector<vtkIdType> face_info(num_fverts);
        for (int fv=0; fv<num_fverts; fv++)
        {
          int v = 0;
          for (int cv=0; cv<num_verts; ++cv)
            if (cell.vertex_ids[cv] == face.vertex_ids[fv])
            { v = cv; break; }

          face_info[fv] = cell_vids[v];
        }

        faces_vids.push_back(static_cast<vtkIdType>(num_fverts));
        for (auto vid : face_info)
          faces_vids.push_back(vid);
      }//for f

      ugrid->InsertNextCell(VTK_POLYHEDRON,
                            static_cast<vtkIdType>(num_verts),
                            cell_vids.data(),
                            static_cast<vtkIdType>(num_faces),
                            faces_vids.data());
    }//polyhedron
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

  for (int g=0; g < ff_uk.num_components; g++)
  {
    ugrid->GetPointData()->AddArray(field_node_array[g]);
    ugrid->GetCellData()->AddArray(field_cell_avg_array[g]);
  }

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  grid_writer->Write();

  //============================================= Parallel summary file
  if (chi_mpi.location_id == 0)
      WritePVTU(base_filename, field_name, static_cast<int>(ff_uk.num_components));
}
