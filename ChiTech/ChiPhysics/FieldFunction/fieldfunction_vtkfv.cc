#include "fieldfunction.h"

#include "ChiMesh/Cell/cell_slab.h"
#include "ChiMesh/Cell/cell_polygon.h"
#include "ChiMesh/Cell/cell_polyhedron.h"
#include "ChiPhysics/chi_physics.h"

#include "ChiMesh/FieldFunctionInterpolation/chi_ffinterpolation.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;
extern ChiPhysics&  chi_physics_handler;

#include <vtkCellType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
//#include <vtkXMLPUnstructuredGridWriter.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include <vtkInformation.h>



//###################################################################
/**Handles the PWLD version of a field function export to VTK.
 *
 * */
void chi_physics::FieldFunction::ExportToVTKFV(const std::string& base_name,
                                               const std::string& field_name)
{
  std::vector<std::vector<double>>    d_nodes;

  auto points = vtkSmartPointer<vtkPoints>::New();

  const auto& field_unknown = unknown_manager.unknowns[ref_variable];

  //============================================= Init grid and material name
  vtkUnstructuredGrid* ugrid;
  vtkIntArray*      matarray;
  vtkIntArray*      pararray;

  vtkDoubleArray*   phiavgarray;

  ugrid    = vtkUnstructuredGrid::New();
  matarray = vtkIntArray::New();
  matarray->SetName("Material");
  pararray = vtkIntArray::New();
  pararray->SetName("Partition");

  phiavgarray = vtkDoubleArray::New();
  phiavgarray->SetName((field_name + std::string("-Avg")).c_str());

  //========================================= Populate nodes
  size_t vcount=0;
  for (const auto& vertex : grid->vertices)
  {
    std::vector<double> d_node;
    d_node.push_back(vertex.x);
    d_node.push_back(vertex.y);
    d_node.push_back(vertex.z);

    d_nodes.push_back(d_node);

    points->InsertPoint(vcount,d_node.data());
    ++vcount;
  }

  //======================================== populate cell mapping
  int num_loc_cells = grid->local_cells.size();
  std::vector<std::pair<uint64_t,uint>> cell_comps_to_map(num_loc_cells);
  std::vector<uint64_t> mapping;
  for (uint64_t lc=0; lc<num_loc_cells; lc++)
    cell_comps_to_map[lc] = std::make_pair(lc,0);

  CreateFVMappingLocal(cell_comps_to_map,
                       mapping);

  //======================================== Populate cell information
  for (const auto& cell : grid->local_cells)
  {
    int lc = cell.local_id;

    int mat_id = cell.material_id;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      auto slab_cell = (chi_mesh::CellSlab*)(&cell);

      std::vector<vtkIdType> cell_info;
      cell_info.push_back(slab_cell->vertex_ids[0]);
      cell_info.push_back(slab_cell->vertex_ids[1]);

      ugrid->
        InsertNextCell(VTK_LINE,2,
                       cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      double phi_value = field_vector_local->operator[](mapping[lc]);
      phiavgarray->InsertNextValue(phi_value);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)(&cell);

      std::vector<vtkIdType> cell_info;

      int num_verts = poly_cell->vertex_ids.size();
      for (int v=0; v<num_verts; v++)
        cell_info.push_back(poly_cell->vertex_ids[v]);

      ugrid->
        InsertNextCell(VTK_POLYGON,num_verts,
                       cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      double phi_value = field_vector_local->operator[](mapping[lc]);
      phiavgarray->InsertNextValue(phi_value);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);

      int num_verts = polyh_cell->vertex_ids.size();
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
        cell_info[v] = polyh_cell->vertex_ids[v];

      vtkSmartPointer<vtkCellArray> faces =
        vtkSmartPointer<vtkCellArray>::New();

      int num_faces = polyh_cell->faces.size();
      for (int f=0; f<num_faces; f++)
      {
        int num_fverts = polyh_cell->faces[f].vertex_ids.size();
        std::vector<vtkIdType> face(num_fverts);
        for (int fv=0; fv<num_fverts; fv++)
          face[fv] = polyh_cell->faces[f].vertex_ids[fv];

        faces->InsertNextCell(num_fverts,face.data());
      }//for f

      ugrid->
        InsertNextCell(VTK_POLYHEDRON,num_verts,
                       cell_info.data(),num_faces,faces->GetPointer());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      double phi_value = field_vector_local->operator[](mapping[lc]);
      phiavgarray->InsertNextValue(phi_value);
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
  vtkXMLUnstructuredGridWriter* grid_writer =
    vtkXMLUnstructuredGridWriter::New();

  ugrid->GetCellData()->AddArray(matarray);
  ugrid->GetCellData()->AddArray(pararray);
  ugrid->GetCellData()->AddArray(phiavgarray);

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  grid_writer->Write();

  //============================================= Parallel summary file
  if (chi_mpi.location_id == 0)
  {
      WritePVTU(base_filename, field_name);
  }
}



//###################################################################
/**Handles the PWLD version of a field function export to VTK.
 *
 * */
void chi_physics::FieldFunction::ExportToVTKFVG(const std::string& base_name,
                                                const std::string& field_name)
{
  chi_mesh::FieldFunctionInterpolation ff_interpol;

  std::vector<std::vector<double>>    d_nodes;

  auto points = vtkSmartPointer<vtkPoints>::New();

  const auto& field_unknown = unknown_manager.unknowns[ref_variable];

  //============================================= Init grid and material name
  vtkUnstructuredGrid* ugrid;
  vtkIntArray*      matarray;
  vtkIntArray*      pararray;

  std::vector<vtkDoubleArray*>   phiavgarray(field_unknown.num_components);

  ugrid    = vtkUnstructuredGrid::New();
  matarray = vtkIntArray::New();
  matarray->SetName("Material");
  pararray = vtkIntArray::New();
  pararray->SetName("Partition");

  for (int g=0; g < field_unknown.num_components; g++)
  {
    phiavgarray[g] = vtkDoubleArray::New();
    phiavgarray[g]->SetName((field_name +
                             std::string("_g") +
                             std::to_string(g) +
                             std::string("-Avg")).c_str());
  }


  //========================================= Populate nodes
  size_t vcount=0;
  for (const auto& vertex : grid->vertices)
  {
    std::vector<double> d_node;
    d_node.push_back(vertex.x);
    d_node.push_back(vertex.y);
    d_node.push_back(vertex.z);

    d_nodes.push_back(d_node);

    points->InsertPoint(vcount,d_node.data());
    ++vcount;
  }

  const auto& ff_uk = this->unknown_manager.unknowns[ref_variable];

  //======================================== populate cell mapping
  int num_loc_cells = grid->local_cells.size();
  std::vector<std::pair<uint64_t,uint>> cell_comps_to_map(num_loc_cells);
  std::vector<uint64_t> mapping;
  for (uint64_t lc=0; lc<num_loc_cells; lc++)
    for (unsigned int g=0; g < ff_uk.num_components; g++)
      cell_comps_to_map[lc] = std::make_pair(lc,g);

  CreateFVMappingLocal(cell_comps_to_map,
                       mapping);


  //======================================== Populate cell information
  int counter=-1;
  for (const auto& cell : grid->local_cells)
  {
    int lc = cell.local_id;

    int mat_id = cell.material_id;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      auto slab_cell = (chi_mesh::CellSlab*)(&cell);

      std::vector<vtkIdType> cell_info;
      cell_info.push_back(slab_cell->vertex_ids[0]);
      cell_info.push_back(slab_cell->vertex_ids[1]);

      ugrid->
        InsertNextCell(VTK_LINE,2,
                       cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      for (int g=0; g < ff_uk.num_components; g++)
      {
        ++counter;
        double phi_value = field_vector_local->operator[](mapping[counter]);
        phiavgarray[g]->InsertNextValue(phi_value);
      }//for g

    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)(&cell);

      std::vector<vtkIdType> cell_info;

      int num_verts = poly_cell->vertex_ids.size();
      for (int v=0; v<num_verts; v++)
        cell_info.push_back(poly_cell->vertex_ids[v]);

      ugrid->
        InsertNextCell(VTK_POLYGON,num_verts,
                       cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      for (int g=0; g < ff_uk.num_components; g++)
      {
        ++counter;
        double phi_value = field_vector_local->operator[](mapping[counter]);
        phiavgarray[g]->InsertNextValue(phi_value);
      }//for g
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);

      int num_verts = polyh_cell->vertex_ids.size();
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
        cell_info[v] = polyh_cell->vertex_ids[v];

      vtkSmartPointer<vtkCellArray> faces =
        vtkSmartPointer<vtkCellArray>::New();

      int num_faces = polyh_cell->faces.size();
      for (int f=0; f<num_faces; f++)
      {
        int num_fverts = polyh_cell->faces[f].vertex_ids.size();
        std::vector<vtkIdType> face(num_fverts);
        for (int fv=0; fv<num_fverts; fv++)
          face[fv] = polyh_cell->faces[f].vertex_ids[fv];

        faces->InsertNextCell(num_fverts,face.data());
      }//for f

      ugrid->
        InsertNextCell(VTK_POLYHEDRON,num_verts,
                       cell_info.data(),num_faces,faces->GetPointer());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      for (int g=0; g < ff_uk.num_components; g++)
      {
        ++counter;
        double phi_value = field_vector_local->operator[](mapping[counter]);
        phiavgarray[g]->InsertNextValue(phi_value);
      }//for g
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
  vtkXMLUnstructuredGridWriter* grid_writer =
    vtkXMLUnstructuredGridWriter::New();

  ugrid->GetCellData()->AddArray(matarray);
  ugrid->GetCellData()->AddArray(pararray);
  for (int g=0; g < ff_uk.num_components; g++)
    ugrid->GetCellData()->AddArray(phiavgarray[g]);

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  /*It seems that cluster systems throw an error when the pvtu file
   * also tries to write to the serial file.*/
//  if (chi_mpi.location_id != 0)
    grid_writer->Write();

  //============================================= Parallel summary file
  if (chi_mpi.location_id == 0)
  {
      WritePVTU(base_filename, field_name, (int)ff_uk.num_components);
  }
}
