#include "fieldfunction.h"

#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>
#include <ChiPhysics/chi_physics.h>

#include <PiecewiseLinear/pwl.h>
#include <PiecewiseLinear/CellViews/pwl_slab.h>
#include <PiecewiseLinear/CellViews/pwl_polygon.h>
#include <PiecewiseLinear/CellViews/pwl_polyhedron.h>

#include <ChiMesh/FieldFunctionInterpolation/chi_ffinterpolation.h>

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog chi_log;
extern ChiMPI chi_mpi;
extern ChiPhysics chi_physics_handler;

#include <vtkCellType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>

#include <vtkInformation.h>

//###################################################################
/**Handles the PWLD version of a field function export to VTK.
 *
 * */
void chi_physics::FieldFunction::ExportToVTKPWLD(const std::string& base_name,
                                                 const std::string& field_name)
{
  auto pwl_sdm = (SpatialDiscretization_PWL*)spatial_discretization;

  chi_mesh::FieldFunctionInterpolation ff_interpol;
  ff_interpol.grid_view = grid;

  std::vector<std::vector<double>>    d_nodes;

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  //============================================= Init grid and material name
  vtkUnstructuredGrid* ugrid;
  vtkIntArray*      matarray;
  vtkIntArray*      pararray;
  vtkDoubleArray*   phiarray;
  vtkDoubleArray*   phiavgarray;

  ugrid    = vtkUnstructuredGrid::New();
  matarray = vtkIntArray::New();
  matarray->SetName("Material");
  pararray = vtkIntArray::New();
  pararray->SetName("Partition");
  phiarray = vtkDoubleArray::New();
  phiarray->SetName(field_name.c_str());
  phiavgarray = vtkDoubleArray::New();
  phiavgarray->SetName((field_name + std::string("-Avg")).c_str());

  //======================================== Populate cell information
  int nc=0;
  for (const auto& cell : grid->local_cells)
  {
    int cell_local_id = cell.local_id;

    int mat_id = cell.material_id;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      auto slab_cell = (chi_mesh::CellSlab*)(&cell);

      int num_verts = 2;
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        int vgi = slab_cell->vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi]->x;
        d_node[1] = grid->vertices[vgi]->y;
        d_node[2] = grid->vertices[vgi]->z;


        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      ugrid->
        InsertNextCell(VTK_LINE,2,
                       cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      //============= Create dof mapping
      std::vector<int> mapping;
      std::vector<int> dofs_to_map(num_verts);
      std::vector<int> cell_to_map(num_verts, cell_local_id);
      for (int v=0; v<num_verts; v++)
        dofs_to_map[v] = v;

      ff_interpol.CreatePWLDMapping(
        num_components,
        num_sets,
        ref_component,
        ref_set,
        dofs_to_map,
        cell_to_map,
        spatial_discretization->cell_dfem_block_address,
        &mapping);

      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        double dof_value = field_vector_local->operator[](mapping[v]);
        cell_avg_value+= dof_value;
        phiarray->InsertNextValue(dof_value);
      }
      phiavgarray->InsertNextValue(cell_avg_value/num_verts);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)(&cell);

      int num_verts = poly_cell->vertex_ids.size();
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        int vgi = poly_cell->vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi]->x;
        d_node[1] = grid->vertices[vgi]->y;
        d_node[2] = grid->vertices[vgi]->z;

        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      ugrid->
        InsertNextCell(VTK_POLYGON,num_verts,
                       cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      //============= Create dof mapping
      std::vector<int> mapping;
      std::vector<int> dofs_to_map(num_verts);
      std::vector<int> cell_to_map(num_verts, cell_local_id);
      for (int v=0; v<num_verts; v++)
        dofs_to_map[v] = v;

      ff_interpol.CreatePWLDMapping(
        num_components,
        num_sets,
        ref_component,
        ref_set,
        dofs_to_map,
        cell_to_map,
        spatial_discretization->cell_dfem_block_address,
        &mapping);

      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        double dof_value = field_vector_local->operator[](mapping[v]);
        cell_avg_value+= dof_value;
        phiarray->InsertNextValue(dof_value);
      }
      phiavgarray->InsertNextValue(cell_avg_value/num_verts);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);
      auto cell_fe_view = (PolyhedronFEView*)pwl_sdm->MapFeViewL(cell.local_id);

      int num_verts = polyh_cell->vertex_ids.size();
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        int vgi = polyh_cell->vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi]->x;
        d_node[1] = grid->vertices[vgi]->y;
        d_node[2] = grid->vertices[vgi]->z;

        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      vtkSmartPointer<vtkCellArray> faces =
        vtkSmartPointer<vtkCellArray>::New();

      int num_faces = polyh_cell->faces.size();
      for (int f=0; f<num_faces; f++)
      {
        int num_fverts = polyh_cell->faces[f].vertex_ids.size();
        std::vector<vtkIdType> face(num_fverts);
        for (int fv=0; fv<num_fverts; fv++)
        {
          int v = cell_fe_view->face_dof_mappings[f][fv];
          face[fv] = cell_info[v];
        }


        faces->InsertNextCell(num_fverts,face.data());
      }//for f

      ugrid->
        InsertNextCell(VTK_POLYHEDRON,num_verts,
                       cell_info.data(),num_faces,faces->GetPointer());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      //============= Create dof mapping
      std::vector<int> mapping;
      std::vector<int> dofs_to_map(num_verts);
      std::vector<int> cell_to_map(num_verts, cell_local_id);
      for (int v=0; v<num_verts; v++)
        dofs_to_map[v] = v;

      ff_interpol.CreatePWLDMapping(
        num_components,
        num_sets,
        ref_component,
        ref_set,
        dofs_to_map,
        cell_to_map,
        spatial_discretization->cell_dfem_block_address,
        &mapping);

      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        double dof_value = field_vector_local->operator[](mapping[v]);
        cell_avg_value+= dof_value;
        phiarray->InsertNextValue(dof_value);
      }
      phiavgarray->InsertNextValue(cell_avg_value/num_verts);

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
  ugrid->GetPointData()->AddArray(phiarray);
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
/**Handles the PWLD version of a field function export to VTK with all groups.
 *
 * */
void chi_physics::FieldFunction::ExportToVTKPWLDG(const std::string& base_name,
                                                  const std::string& field_name)
{
  SpatialDiscretization_PWL* pwl_sdm =
    (SpatialDiscretization_PWL*)spatial_discretization;

  chi_mesh::FieldFunctionInterpolation ff_interpol;
  ff_interpol.grid_view = grid;

  std::vector<std::vector<double>>    d_nodes;

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  //============================================= Init grid and material name
  vtkUnstructuredGrid* ugrid;
  vtkIntArray*      matarray;
  vtkIntArray*      pararray;
  std::vector<vtkDoubleArray*>   phiarray(num_components);
  std::vector<vtkDoubleArray*>   phiavgarray(num_components);

  ugrid    = vtkUnstructuredGrid::New();
  matarray = vtkIntArray::New();
  matarray->SetName("Material");
  pararray = vtkIntArray::New();
  pararray->SetName("Partition");

  for (int g=0; g < num_components; g++)
  {
    char group_text[100];
    sprintf(group_text,"%03d",g);
    phiarray[g]    = vtkDoubleArray::New();
    phiavgarray[g] = vtkDoubleArray::New();

    phiarray[g]   ->SetName((field_name +
                             std::string("_g") +
                             std::string(group_text)).c_str());
    phiavgarray[g]->SetName((field_name +
                             std::string("_g") +
                             std::string(group_text) +
                             std::string("_avg")).c_str());
  }


  //======================================== Populate cell information
  int nc=0;
  for (const auto& cell : grid->local_cells)
  {
    int cell_local_id = cell.local_id;

    int mat_id = cell.material_id;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      auto slab_cell = (chi_mesh::CellSlab*)(&cell);

      int num_verts = 2;
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        int vgi = slab_cell->vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi]->x;
        d_node[1] = grid->vertices[vgi]->y;
        d_node[2] = grid->vertices[vgi]->z;


        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      ugrid->
        InsertNextCell(VTK_LINE,2,
                       cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      //============= Create dof mapping
      std::vector<int> mapping;
      std::vector<int> dofs_to_map(num_verts);
      std::vector<int> cell_to_map(num_verts, cell_local_id);
      for (int v=0; v<num_verts; v++)
        dofs_to_map[v] = v;

      ff_interpol.CreatePWLDMapping(
        num_components,
        num_sets,
        ref_component,
        ref_set,
        dofs_to_map,
        cell_to_map,
        spatial_discretization->cell_dfem_block_address,
        &mapping);

      for (int g=0; g < num_components; g++)
      {
        double cell_avg_value = 0.0;
        for (int v=0; v<num_verts; v++)
        {
          double dof_value = field_vector_local->operator[](mapping[v]+g);
          cell_avg_value+= dof_value;
          phiarray[g]->InsertNextValue(dof_value);
        }
        phiavgarray[g]->InsertNextValue(cell_avg_value/num_verts);
      }//for g

    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)(&cell);

      int num_verts = poly_cell->vertex_ids.size();
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        int vgi = poly_cell->vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi]->x;
        d_node[1] = grid->vertices[vgi]->y;
        d_node[2] = grid->vertices[vgi]->z;

        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      ugrid->
        InsertNextCell(VTK_POLYGON,num_verts,
                       cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      //============= Create dof mapping
      std::vector<int> mapping;
      std::vector<int> dofs_to_map(num_verts);
      std::vector<int> cell_to_map(num_verts, cell_local_id);
      for (int v=0; v<num_verts; v++)
        dofs_to_map[v] = v;

      ff_interpol.CreatePWLDMapping(
        num_components,
        num_sets,
        ref_component,
        ref_set,
        dofs_to_map,
        cell_to_map,
        spatial_discretization->cell_dfem_block_address,
        &mapping);

      for (int g=0; g < num_components; g++)
      {
        double cell_avg_value = 0.0;
        for (int v=0; v<num_verts; v++)
        {
          double dof_value = field_vector_local->operator[](mapping[v]+g);
          cell_avg_value+= dof_value;
          phiarray[g]->InsertNextValue(dof_value);
        }
        phiavgarray[g]->InsertNextValue(cell_avg_value/num_verts);
      }//for g
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);
      auto cell_fe_view = (PolyhedronFEView*)pwl_sdm->MapFeViewL(cell.local_id);

      int num_verts = polyh_cell->vertex_ids.size();
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        int vgi = polyh_cell->vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi]->x;
        d_node[1] = grid->vertices[vgi]->y;
        d_node[2] = grid->vertices[vgi]->z;

        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      vtkSmartPointer<vtkCellArray> faces =
        vtkSmartPointer<vtkCellArray>::New();

      int num_faces = polyh_cell->faces.size();
      for (int f=0; f<num_faces; f++)
      {
        int num_fverts = polyh_cell->faces[f].vertex_ids.size();
        std::vector<vtkIdType> face(num_fverts);
        for (int fv=0; fv<num_fverts; fv++)
        {
          int v = cell_fe_view->face_dof_mappings[f][fv];
          face[fv] = cell_info[v];
        }


        faces->InsertNextCell(num_fverts,face.data());
      }//for f

      ugrid->
        InsertNextCell(VTK_POLYHEDRON,num_verts,
                       cell_info.data(),num_faces,faces->GetPointer());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      //============= Create dof mapping
      std::vector<int> mapping;
      std::vector<int> dofs_to_map(num_verts);
      std::vector<int> cell_to_map(num_verts, cell_local_id);
      for (int v=0; v<num_verts; v++)
        dofs_to_map[v] = v;

      ff_interpol.CreatePWLDMapping(
        num_components,
        num_sets,
        ref_component,
        ref_set,
        dofs_to_map,
        cell_to_map,
        spatial_discretization->cell_dfem_block_address,
        &mapping);

      for (int g=0; g < num_components; g++)
      {
        double cell_avg_value = 0.0;
        for (int v=0; v<num_verts; v++)
        {
          double dof_value = field_vector_local->operator[](mapping[v]+g);
          cell_avg_value+= dof_value;
          phiarray[g]->InsertNextValue(dof_value);
        }
        phiavgarray[g]->InsertNextValue(cell_avg_value/num_verts);
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

  for (int g=0; g < num_components; g++)
  {
    ugrid->GetPointData()->AddArray(phiarray[g]);
    ugrid->GetCellData()->AddArray(phiavgarray[g]);
  }


  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  /*It seems that cluster systems throw an error when the pvtu file
   * also tries to write to the serial file.
  if (chi_mpi.location_id != 0)
   */
    grid_writer->Write();


  //============================================= Parallel summary file
  if (chi_mpi.location_id == 0)
  {
      WritePVTU(base_filename, field_name, num_components);
  }
}