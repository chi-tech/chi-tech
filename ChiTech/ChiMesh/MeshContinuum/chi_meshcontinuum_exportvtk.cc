#include "chi_meshcontinuum.h"

#include "ChiMesh/Cell/cell.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include <vtkCellType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include <vtkInformation.h>

//###################################################################
/**Exports just the mesh to VTK format.*/
void chi_mesh::MeshContinuum::ExportCellsToVTK(const char* baseName) const
{
  chi_log.Log() << "Exporting mesh to VTK. " << local_cells.size();
  std::vector<std::vector<double>> d_nodes;

  auto points = vtkSmartPointer<vtkPoints>::New();

  //============================================= Init grid and material name
  auto ugrid    = vtkUnstructuredGrid::New();
  auto matarray = vtkIntArray::New();
  auto pararray = vtkIntArray::New();

  matarray->SetName("Material");
  pararray->SetName("Partition");

  auto grid = this;

  //======================================== Populate cell information
  size_t cell_count = 0;
  int64_t node_count = 0;
  for (const auto& cell : grid->local_cells)
  {
    int mat_id = cell.material_id;

    const size_t num_verts = cell.vertex_ids.size();
    std::vector<vtkIdType> cell_vids(num_verts);

    for (size_t v=0; v<num_verts; ++v)
    {
      uint64_t vgi = cell.vertex_ids[v]; //vertex global id
      std::vector<double> d_node(3);
      auto vertex = grid->vertices[vgi];
      d_node[0] = vertex.x;
      d_node[1] = vertex.y;
      d_node[2] = vertex.z;

      points->InsertPoint(node_count, d_node.data());
      cell_vids[v] = node_count++;
    }

    matarray->InsertNextValue(mat_id);
    pararray->InsertNextValue(static_cast<int>(cell.partition_id));

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell.Type() == chi_mesh::CellType::SLAB)
      ugrid->InsertNextCell(VTK_LINE, 2, cell_vids.data());

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
      ugrid->InsertNextCell(VTK_POLYGON,
                            static_cast<vtkIdType>(num_verts),
                            cell_vids.data());

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      // Build polyhedron faces
      std::vector<vtkIdType> faces_vids;

      size_t num_faces = cell.faces.size();
      for (auto& face : cell.faces)
      {
        size_t num_fverts = face.vertex_ids.size();
        std::vector<vtkIdType> face_info(num_fverts);
        for (size_t fv=0; fv<num_fverts; fv++)
        {
          size_t v = 0;
          for (size_t cv=0; cv<num_verts; ++cv)
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
    ++cell_count;
  }//for local cells

  ugrid->SetPoints(points);

  //============================================= Construct file name
  std::string base_filename     = std::string(baseName);
  std::string location_filename = base_filename +
                                  std::string("_") +
                                  std::to_string(chi_mpi.location_id) +
                                  std::string(".vtu");

  //============================================= Serial Output each piece
  auto grid_writer = vtkXMLUnstructuredGridWriter::New();

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

    // Cut off path to base_filename
    std::string filename_short = base_filename.substr(base_filename.find_last_of("/\\")+1);

    for (int p=0; p<chi_mpi.process_count; ++p)
    {
      if (is_global_mesh and p!=0) continue;

      ofile << "      <Piece Source=\""
            << filename_short +
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