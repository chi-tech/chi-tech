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
void chi_physics::FieldFunction::ExportToVTKPWLC(const std::string& base_name,
                                                 const std::string& field_name)
{
  SpatialDiscretization_PWL* pwl_sdm =
    (SpatialDiscretization_PWL*)spatial_discretization;

  chi_mesh::FieldFunctionInterpolation ff_interpol;

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

  //======================================== Precreate nodes to map
  std::vector<int> cfem_nodes;

  for (const auto& cell : grid->local_cells)
  {
    int cell_g_ind = cell.cell_global_id;

    int mat_id = cell.material_id;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      auto slab_cell = (chi_mesh::CellSlab*)(&cell);

      int num_verts = 2;
      for (int v=0; v<num_verts; v++)
        cfem_nodes.push_back(slab_cell->vertex_ids[v]);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)(&cell);

      int num_verts = poly_cell->vertex_ids.size();
      for (int v=0; v<num_verts; v++)
        cfem_nodes.push_back(poly_cell->vertex_ids[v]);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);

      int num_verts = polyh_cell->vertex_ids.size();
      for (int v=0; v<num_verts; v++)
        cfem_nodes.push_back(polyh_cell->vertex_ids[v]);
    }//polyhedron
  }//for local cells

  std::vector<int> mapping;
  Vec phi_vec;
  ff_interpol.CreateCFEMMapping(num_components, num_sets, ref_component, ref_set,
                                *field_vector, phi_vec, cfem_nodes, &mapping,
                                spatial_discretization);


  //======================================== Populate cell information
  int nc=0;
  int counter=-1;
  for (const auto& cell : grid->local_cells)
  {
    int cell_g_ind = cell.cell_global_id;

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

      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        counter++;
        int ir = mapping[counter];
        double dof_value = 0.0;
        VecGetValues(phi_vec,1,&ir,&dof_value);;
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

      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        counter++;
        int ir = mapping[counter];
        double dof_value = 0.0;
        VecGetValues(phi_vec,1,&ir,&dof_value);;
        cell_avg_value+= dof_value;
        phiarray->InsertNextValue(dof_value);
      }
      phiavgarray->InsertNextValue(cell_avg_value/num_verts);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);
      auto cell_fe_view = (PolyhedronFEView*)pwl_sdm->MapFeViewL(cell.cell_local_id);

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

      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        counter++;
        int ir = mapping[counter];
        double dof_value = 0.0;
        VecGetValues(phi_vec,1,&ir,&dof_value);;
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
/**Handles the PWLD version of a field function export to VTK.
 *
 * */
void chi_physics::FieldFunction::ExportToVTKPWLCG(const std::string& base_name,
                                                  const std::string& field_name)
{
  SpatialDiscretization_PWL* pwl_sdm =
    (SpatialDiscretization_PWL*)spatial_discretization;

  chi_mesh::FieldFunctionInterpolation ff_interpol;

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

  //======================================== Precreate nodes to map
  std::vector<int> cfem_nodes;

  for (const auto& cell : grid->local_cells)
  {
    int cell_g_ind = cell.cell_global_id;

    int mat_id = cell.material_id;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      auto slab_cell = (chi_mesh::CellSlab*)(&cell);

      int num_verts = 2;
      for (int v=0; v<num_verts; v++)
        for (int g=0; g < num_components; g++)
          cfem_nodes.push_back(slab_cell->vertex_ids[v]);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)(&cell);

      int num_verts = poly_cell->vertex_ids.size();
      for (int v=0; v<num_verts; v++)
        for (int g=0; g < num_components; g++)
          cfem_nodes.push_back(poly_cell->vertex_ids[v]);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);

      int num_verts = polyh_cell->vertex_ids.size();
      for (int v=0; v<num_verts; v++)
        for (int g=0; g < num_components; g++)
          cfem_nodes.push_back(polyh_cell->vertex_ids[v]);
    }//polyhedron
  }//for local cells

  std::vector<int> mapping;
  Vec phi_vec;
  ff_interpol.CreateCFEMMapping(num_components, num_sets, ref_component, ref_set,
                                *field_vector, phi_vec, cfem_nodes, &mapping,
                                spatial_discretization);


  //======================================== Populate cell information
  int nc=0;
  int counter=-1;
  for (const auto& cell : grid->local_cells)
  {
    int cell_g_ind = cell.cell_global_id;

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

      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        counter++;
        int ir = mapping[counter];
        double dof_value = 0.0;
        VecGetValues(phi_vec,1,&ir,&dof_value);;
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

      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        counter++;
        int ir = mapping[counter];
        double dof_value = 0.0;
        VecGetValues(phi_vec,1,&ir,&dof_value);;
        cell_avg_value+= dof_value;
        phiarray->InsertNextValue(dof_value);
      }
      phiavgarray->InsertNextValue(cell_avg_value/num_verts);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);
      auto cell_fe_view = (PolyhedronFEView*)pwl_sdm->MapFeViewL(cell.cell_local_id);

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

      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        counter++;
        int ir = mapping[counter];
        double dof_value = 0.0;
        VecGetValues(phi_vec,1,&ir,&dof_value);;
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

  /*It seems that cluster systems throw an error when the pvtu file
   * also tries to write to the serial file.*/
  if (chi_mpi.location_id != 0)
    grid_writer->Write();

  //============================================= Parallel summary file
  if (chi_mpi.location_id == 0)
  {
      WritePVTU(base_filename, field_name, num_components);
  }
}