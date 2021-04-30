#include "chi_meshcontinuum.h"

#include "ChiMesh/Cell/cell_slab.h"
#include "ChiMesh/Cell/cell_polygon.h"
#include "ChiMesh/Cell/cell_polyhedron.h"
#include <ChiPhysics/chi_physics.h>

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;
extern ChiPhysics&  chi_physics_handler;

#include <vtkCellType.h>
#include <vtkUnstructuredGrid.h>
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
/**On hold*/
void chi_mesh::MeshContinuum::ExportCellsToVTK(const char* baseName)
{
  chi_log.Log() << "Exporting mesh to VTK. " << local_cells.size();
  std::vector<std::vector<double>>    d_nodes;

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  //============================================= Init grid and material name
  auto ugrid    = vtkUnstructuredGrid::New();
  auto matarray = vtkIntArray::New();
  auto pararray = vtkIntArray::New();

  matarray->SetName("Material");
  pararray->SetName("Partition");

  //======================================== Precreate nodes to map
  std::vector<int> cfem_nodes;

  auto grid = this;

  for (const auto& cell : grid->local_cells)
    for (auto vid : cell.vertex_ids)
      cfem_nodes.push_back(vid);

  //======================================== Populate cell information
  int nc=0;
  for (const auto& cell : grid->local_cells)
  {
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
        d_node[0] = grid->vertices[vgi].x;
        d_node[1] = grid->vertices[vgi].y;
        d_node[2] = grid->vertices[vgi].z;


        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      ugrid->InsertNextCell(VTK_LINE, 2, cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)(&cell);

      size_t num_verts = poly_cell->vertex_ids.size();
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        uint64_t vgi = poly_cell->vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi].x;
        d_node[1] = grid->vertices[vgi].y;
        d_node[2] = grid->vertices[vgi].z;

        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      ugrid->InsertNextCell(VTK_POLYGON, num_verts, cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);

      size_t num_verts = polyh_cell->vertex_ids.size();
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        uint64_t vgi = polyh_cell->vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi].x;
        d_node[1] = grid->vertices[vgi].y;
        d_node[2] = grid->vertices[vgi].z;

        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      vtkSmartPointer<vtkCellArray> faces =
        vtkSmartPointer<vtkCellArray>::New();

      size_t num_faces = polyh_cell->faces.size();
      for (int f=0; f<num_faces; f++)
      {
        size_t num_fverts = polyh_cell->faces[f].vertex_ids.size();
        std::vector<vtkIdType> face(num_fverts);
        for (int fv=0; fv<num_fverts; fv++)
        {
          int v = -1;
          for (int vr=0; vr<cell.vertex_ids.size(); ++vr)
            if (polyh_cell->faces[f].vertex_ids[fv] ==
                cell.vertex_ids[vr])
            {
              v = vr; break;
            }
          face[fv] = cell_info[v];
        }


        faces->InsertNextCell(num_fverts,face.data());
      }//for f

      ugrid->
        InsertNextCell(VTK_POLYHEDRON,num_verts,
                       cell_info.data(),num_faces,faces->GetPointer());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);
    }//polyhedron
  }//for local cells

  ugrid->SetPoints(points);

  //============================================= Construct file name
  std::string base_filename     = std::string(baseName);
  std::string location_filename = base_filename +
                                  std::string("_") +
                                  std::to_string(chi_mpi.location_id) +
                                  std::string(".vtu");

  //============================================= Serial Output each piece
  vtkXMLUnstructuredGridWriter* grid_writer =
    vtkXMLUnstructuredGridWriter::New();

  ugrid->GetCellData()->AddArray(matarray);
  ugrid->GetCellData()->AddArray(pararray);

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  grid_writer->Write();

  //============================================= Parallel summary file
  if (chi_mpi.location_id == 0)
  {
    std::string summary_file_name = base_filename + std::string(".pvtu");
    std::ofstream ofile;
    ofile.open(summary_file_name);

    ofile << "<?xml version=\"1.0\"?>" << std::endl;
    ofile << "<!--" << std::endl;
    ofile << "#Unstructured Mesh" << std::endl;
    ofile << "-->" << std::endl;
    ofile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
          << "byte_order=\"LittleEndian\">" << std::endl;
    ofile << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
    ofile << "    <PCellData Scalars=\"scalars\">" << std::endl;
    ofile << "      <PDataArray type=\"Int32\" Name=\"Material\" "
          << " format=\"ascii\"/>" << std::endl;
    ofile << "      <PDataArray type=\"Int32\" Name=\"Partition\""
          << " format=\"ascii\"/>" << std::endl;

    ofile << "    </PCellData>" << std::endl;
    ofile << "    <PPoints>" << std::endl;
    ofile << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>" << std::endl;
    ofile << "    </PPoints>" << std::endl;

    bool is_global_mesh =
      chi_mesh::GetCurrentHandler()->volume_mesher->options.mesh_global;

    for (int p=0; p<chi_mpi.process_count; p++)
    {
      if (is_global_mesh and p!=0) continue;

      ofile << "      <Piece Source=\""
            << base_filename +
               std::string("_") +
               std::to_string(p) +
               std::string(".vtu")
            << "\"/>" << std::endl;
    }

    ofile << "  </PUnstructuredGrid>" << std::endl;
    ofile << "</VTKFile>" << std::endl;

    ofile.close();
  }

  chi_log.Log() << "Done exporting mesh to VTK.";
}