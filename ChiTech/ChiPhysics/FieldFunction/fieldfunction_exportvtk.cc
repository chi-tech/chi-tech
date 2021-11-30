#include "fieldfunction.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include <vtkCellType.h>

#include <fstream>

//###################################################################
/**Exports a field function to VTK format.
 *
 * */
void chi_physics::FieldFunction::
  ExportToVTKComponentOnly(const std::string& base_name,
                           const std::string& field_name)
{
  chi_log.Log(LOG_0)
    << "Exporting field function " << text_name
    << " to files with base name " << base_name;

  typedef chi_math::SpatialDiscretizationType SDMType;
  auto& field_sdm_type = spatial_discretization->type;

  if (field_sdm_type == SDMType::FINITE_VOLUME)
    ExportToVTKFV(base_name,field_name);
  if (field_sdm_type == SDMType::PIECEWISE_LINEAR_CONTINUOUS)
    ExportToVTKPWLC(base_name,field_name);
  if (field_sdm_type == SDMType::PIECEWISE_LINEAR_DISCONTINUOUS)
    ExportToVTKPWLD(base_name,field_name);

}

//###################################################################
/**Exports a field function to VTK format but exports all of the
 * available groups.
 *
 * */
void chi_physics::FieldFunction::ExportToVTK(const std::string& base_name,
                                             const std::string& field_name)
{
  chi_log.Log(LOG_0)
    << "Exporting field function " << text_name
    << " to files with base name " << base_name
    << " to field name " << field_name;

  typedef chi_math::SpatialDiscretizationType SDMType;
  auto& field_sdm_type = spatial_discretization->type;

  if (field_sdm_type == SDMType::FINITE_VOLUME)
    ExportToVTKFV(base_name, field_name, true);
  if (field_sdm_type == SDMType::PIECEWISE_LINEAR_CONTINUOUS)
    ExportToVTKPWLC(base_name, field_name, true);
  if (field_sdm_type == SDMType::PIECEWISE_LINEAR_DISCONTINUOUS)
    ExportToVTKPWLD(base_name, field_name, true);

}


//###################################################################
/** Writes the VTK "Assembly file" for multiple vtu files.

\param base_filename Base name for all vtu file. .pvtu will get appended to it.
\param field_name Name of the field to be exported
\param num_components Optional. Defaults to 0. If greater than zero then
                        exports components.
\author Jason*/
void chi_physics::FieldFunction::WritePVTU(const std::string& base_filename,
                                           const std::string& field_name,
                                           const std::vector<std::string>& component_names)
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
  ofile << "    <PPointData Scalars=\"scalars\">" << std::endl;

  const size_t num_components = component_names.size();

  for (size_t c=0; c < num_components; c++)
  {
    ofile << "      <PDataArray type=\"Float64\" Name=\""
          << component_names[c]
          << "\" format=\"ascii\"/>" << std::endl;
  }

  ofile << "    </PPointData>" << std::endl;
  ofile << "    <PCellData Scalars=\"scalars\">" << std::endl;
  ofile << "      <PDataArray type=\"Int32\" Name=\"Material\" "
        << " format=\"ascii\"/>" << std::endl;
  ofile << "      <PDataArray type=\"UInt32\" Name=\"Partition\""
        << " format=\"ascii\"/>" << std::endl;

  for (size_t c=0; c < num_components; c++)
  {
    ofile << "      <PDataArray type=\"Float64\" Name=\""
          << component_names[c] + std::string("-avg")
          << "\" format=\"ascii\"/>" << std::endl;
  }

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


//###################################################################
/**Uploads just the geometry portion of a cell to VTK.*/
void chi_physics::FieldFunction::
  UploadCellGeometry(const chi_mesh::Cell &cell,
                     int64_t& node_counter,
                     vtkNew<vtkPoints>& points,
                     vtkNew<vtkUnstructuredGrid> &ugrid)
{
  size_t num_verts = cell.vertex_ids.size();

  std::vector<vtkIdType> cell_vids(num_verts);
  for (size_t v=0; v<num_verts; v++)
  {
    uint64_t vgi = cell.vertex_ids[v];
    std::vector<double> d_node(3);
    d_node[0] = grid->vertices[vgi].x;
    d_node[1] = grid->vertices[vgi].y;
    d_node[2] = grid->vertices[vgi].z;

    points->InsertPoint(node_counter,d_node.data());
    cell_vids[v] = node_counter++;
  }

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
}
