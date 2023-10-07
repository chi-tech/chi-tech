#include "SplitFileMeshGenerator.h"

#include "data_types/byte_array.h"
#include "utils/chi_utils.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/VolumeMesher/chi_volumemesher.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
#include "utils/chi_timer.h"

#include "ChiObjectFactory.h"

#include <filesystem>

#define scint static_cast<int>

namespace chi_mesh
{

RegisterChiObject(chi_mesh, SplitFileMeshGenerator);

// ##################################################################
chi::InputParameters SplitFileMeshGenerator::GetInputParameters()
{
  chi::InputParameters params = MeshGenerator::GetInputParameters();

  params.SetGeneralDescription(
    "Generates the mesh only on location 0, thereafter partitions the mesh"
    " but instead of broadcasting the mesh to other locations it creates binary"
    " mesh files for each location.");
  params.SetDocGroup("doc_MeshGenerators");

  params.AddOptionalParameter("num_parts",
                              0,
                              "The number of parts to generate. If zero will "
                              "default to the number of MPI processes. Is "
                              "ignored if the number of MPI processes > 1.");

  params.AddOptionalParameter(
    "split_mesh_dir_path",
    "SplitMesh",
    "Path of the directory to be created for containing the split meshes.");

  params.AddOptionalParameter("split_file_prefix",
                              "split_mesh",
                              "Prefix to use for all split mesh files");

  params.AddOptionalParameter(
    "read_only",
    false,
    "Controls whether the split mesh is recreated or just read.");

  return params;
}

// ##################################################################
SplitFileMeshGenerator::SplitFileMeshGenerator(
  const chi::InputParameters& params)
  : MeshGenerator(params),
    num_parts_(params.GetParamValue<int>("num_parts")),
    split_mesh_dir_path_(
      params.GetParamValue<std::string>("split_mesh_dir_path")),
    split_file_prefix_(params.GetParamValue<std::string>("split_file_prefix")),
    read_only_(params.GetParamValue<bool>("read_only"))
{
}

// ##################################################################
void SplitFileMeshGenerator::Execute()
{
  const int num_mpi = Chi::mpi.process_count;
  const int num_parts = num_mpi == 1 ? num_parts_ : num_mpi;

  if (Chi::mpi.location_id == 0 and (not read_only_))
  {
    //======================================== Execute all input generators
    // Note these could be empty
    std::unique_ptr<UnpartitionedMesh> current_umesh = nullptr;
    for (auto mesh_generator_ptr : inputs_)
    {
      auto new_umesh =
        mesh_generator_ptr->GenerateUnpartitionedMesh(std::move(current_umesh));
      current_umesh = std::move(new_umesh);
    }

    //======================================== Generate final umesh
    current_umesh = GenerateUnpartitionedMesh(std::move(current_umesh));

    Chi::log.Log() << "Writing split-mesh with " << num_parts << " parts";
    const auto cell_pids = PartitionMesh(*current_umesh, num_parts);
    WriteSplitMesh(cell_pids, *current_umesh, num_parts);
    Chi::log.Log() << "Split-mesh with " << num_parts
                   << " parts successfully created";
  } // if home location

  // Other locations wait here for files to be written
  Chi::mpi.Barrier();

  if (Chi::mpi.process_count == num_parts)
  {
    Chi::log.Log() << "Reading split-mesh";
    auto mesh_info = ReadSplitMesh();

    auto grid_ptr = SetupLocalMesh(mesh_info);

    auto new_mesher =
      std::make_shared<chi_mesh::VolumeMesher>(VolumeMesherType::UNPARTITIONED);
    new_mesher->SetContinuum(grid_ptr);

    if (Chi::current_mesh_handler < 0) chi_mesh::PushNewHandlerAndGetIndex();

    auto& cur_hndlr = chi_mesh::GetCurrentHandler();
    cur_hndlr.SetVolumeMesher(new_mesher);
    Chi::log.Log() << "Done reading split-mesh files";
  }
  else
  {
    Chi::log.Log0Warning()
      << "After creating a split-mesh with mpi-processes < "
         "num_parts the program will now auto terminate. This is not an error "
         "and is the default behavior for the SplitFileMeshGenerator.";
    Chi::Exit(EXIT_SUCCESS);
  }

  Chi::mpi.Barrier();
}

// ##################################################################
void SplitFileMeshGenerator::WriteSplitMesh(
  const std::vector<int64_t>& cell_pids,
  const UnpartitionedMesh& umesh,
  int num_parts)
{
  const std::filesystem::path dir_path =
    std::filesystem::absolute(split_mesh_dir_path_);

  const auto parent_path = dir_path.parent_path();
  ChiInvalidArgumentIf(not std::filesystem::exists(parent_path),
                       "Parent path " + parent_path.string() +
                         " does not exist");

  bool root_dir_created = true;
  if (not std::filesystem::exists(dir_path))
    root_dir_created = std::filesystem::create_directories(dir_path);

  ChiLogicalErrorIf(not root_dir_created,
                    "Failed to create directory " + dir_path.string());

  const auto& vertex_subs = umesh.GetVertextCellSubscriptions();
  const auto& raw_cells = umesh.GetRawCells();
  const auto& raw_vertices = umesh.GetVertices();

  for (int pid = 0; pid < num_parts; ++pid)
  {
    const std::filesystem::path file_path = dir_path.string() + "/" +
                                            split_file_prefix_ + "_" +
                                            std::to_string(pid) + ".cmesh";

    std::ofstream ofile(file_path.string(),
                        std::ios_base::binary | std::ios_base::out);

    ChiLogicalErrorIf(not ofile.is_open(),
                      "Failed to open " + file_path.string());

    // Appropriate cells and vertices to the current part being writting
    std::vector<CellPIDGID> cells_needed;
    std::set<uint64_t> vertices_needed;
    {
      uint64_t cell_global_id = 0;
      for (const auto& raw_cell : raw_cells)
      {
        const auto cell_pid = cell_pids[cell_global_id];
        if (CellHasLocalScope(
              pid, *raw_cell, cell_global_id, vertex_subs, cell_pids))
        {
          cells_needed.emplace_back(cell_pid, cell_global_id);

          for (uint64_t vid : raw_cell->vertex_ids)
            vertices_needed.insert(vid);
        }
        ++cell_global_id;
      } // for raw cell
    }

    //================================================ Write mesh attributes
    //                                                 and general info
    const auto& mesh_options = umesh.GetMeshOptions();

    chi::WriteBinaryValue(ofile, num_parts); // int

    chi::WriteBinaryValue(ofile, scint(umesh.GetMeshAttributes())); // int
    chi::WriteBinaryValue(ofile, mesh_options.ortho_Nx);            // size_t
    chi::WriteBinaryValue(ofile, mesh_options.ortho_Ny);            // size_t
    chi::WriteBinaryValue(ofile, mesh_options.ortho_Nz);            // size_t

    chi::WriteBinaryValue(ofile, raw_vertices.size()); // size_t

    //================================================ Write the boundary map
    const auto& bndry_map = mesh_options.boundary_id_map;
    chi::WriteBinaryValue(ofile, bndry_map.size()); // size_t
    for (const auto& [bid, bname] : bndry_map)
    {
      chi::WriteBinaryValue(ofile, bid); // uint64_t
      const size_t num_chars = bname.size();
      chi::WriteBinaryValue(ofile, num_chars);     // size_t
      ofile.write(bname.data(), scint(num_chars)); // characters
    }

    //================================================ Write how many cells
    //                                                 and vertices in file
    chi::WriteBinaryValue(ofile, cells_needed.size());    // size_t
    chi::WriteBinaryValue(ofile, vertices_needed.size()); // size_t

    //================================================ Write cells
    for (const auto& [cell_pid, cell_global_id] : cells_needed)
    {
      const auto& cell = *raw_cells.at(cell_global_id);
      chi_data_types::ByteArray serial_data;
      serial_data.Write(cell_pid);       // int
      serial_data.Write(cell_global_id); // uint64_t
      SerializeCell(cell, serial_data);
      ofile.write((char*)serial_data.Data().data(), scint(serial_data.Size()));
    }

    //================================================ Write vertices
    for (const uint64_t vid : vertices_needed)
    {
      chi_data_types::ByteArray serial_data;
      serial_data.Write(vid); // uint64_t
      serial_data.Write(raw_vertices[vid]);
      ofile.write((char*)serial_data.Data().data(), scint(serial_data.Size()));
    }

    ofile.close();
  } // for p
}

// ##################################################################
void SplitFileMeshGenerator::SerializeCell(
  const UnpartitionedMesh::LightWeightCell& cell,
  chi_data_types::ByteArray& serial_buffer)
{
  serial_buffer.Write(cell.type);
  serial_buffer.Write(cell.sub_type);
  serial_buffer.Write(cell.centroid);
  serial_buffer.Write(cell.material_id);
  serial_buffer.Write(cell.vertex_ids.size());
  for (uint64_t vid : cell.vertex_ids)
    serial_buffer.Write(vid);
  serial_buffer.Write(cell.faces.size());
  for (const auto& face : cell.faces)
  {
    serial_buffer.Write(face.vertex_ids.size());
    for (uint64_t vid : face.vertex_ids)
      serial_buffer.Write(vid);
    serial_buffer.Write(face.has_neighbor);
    serial_buffer.Write(face.neighbor);
  }
}

SplitFileMeshGenerator::SplitMeshInfo SplitFileMeshGenerator::ReadSplitMesh()
{
  const int pid = Chi::mpi.location_id;
  const std::filesystem::path dir_path =
    std::filesystem::absolute(split_mesh_dir_path_);
  const std::filesystem::path file_path = dir_path.string() + "/" +
                                          split_file_prefix_ + "_" +
                                          std::to_string(pid) + ".cmesh";

  SplitMeshInfo info_block;
  auto& cells = info_block.cells_;
  auto& vertices = info_block.vertices_;
  std::ifstream ifile(file_path, std::ios_base::binary | std::ios_base::in);

  ChiLogicalErrorIf(not ifile.is_open(),
                    "Failed to open " + file_path.string());

  //================================================== Read mesh attributes
  //                                                   and general info
  const size_t file_num_parts = chi::ReadBinaryValue<int>(ifile);

  ChiLogicalErrorIf(Chi::mpi.process_count != file_num_parts,
                    "Split mesh files with prefix \"" + split_file_prefix_ +
                      "\" has been created with " +
                      std::to_string(file_num_parts) +
                      " parts but is now being read with " +
                      std::to_string(Chi::mpi.process_count) + " processes.");

  info_block.mesh_attributes_ = chi::ReadBinaryValue<int>(ifile);
  info_block.ortho_Nx_ = chi::ReadBinaryValue<size_t>(ifile);
  info_block.ortho_Ny_ = chi::ReadBinaryValue<size_t>(ifile);
  info_block.ortho_Nz_ = chi::ReadBinaryValue<size_t>(ifile);

  info_block.num_global_vertices_ = chi::ReadBinaryValue<size_t>(ifile);

  //================================================== Read boundary map
  const size_t num_bndries = chi::ReadBinaryValue<size_t>(ifile);
  for (size_t b = 0; b < num_bndries; ++b)
  {
    const uint64_t bid = chi::ReadBinaryValue<uint64_t>(ifile);
    const size_t num_chars = chi::ReadBinaryValue<size_t>(ifile);
    std::string bname(num_chars, ' ');
    ifile.read(bname.data(), static_cast<int>(num_chars));

    info_block.boundary_id_map_.insert(std::make_pair(bid, bname));
  }

  //================================================ Write how many cells
  //                                                 and vertices in file
  const size_t num_cells = chi::ReadBinaryValue<size_t>(ifile);
  const size_t num_vertices = chi::ReadBinaryValue<size_t>(ifile);

  //================================================== Read the cells
  for (size_t c = 0; c < num_cells; ++c)
  {
    const int cell_pid = chi::ReadBinaryValue<int>(ifile);
    const uint64_t cell_gid = chi::ReadBinaryValue<uint64_t>(ifile);
    const CellType cell_type = chi::ReadBinaryValue<CellType>(ifile);
    const CellType cell_sub_type = chi::ReadBinaryValue<CellType>(ifile);

    UnpartitionedMesh::LightWeightCell new_cell(cell_type, cell_sub_type);

    new_cell.centroid = chi::ReadBinaryValue<chi_mesh::Vector3>(ifile);
    new_cell.material_id = chi::ReadBinaryValue<int>(ifile);

    const size_t num_vids = chi::ReadBinaryValue<size_t>(ifile);
    for (size_t v = 0; v < num_vids; ++v)
      new_cell.vertex_ids.push_back(chi::ReadBinaryValue<uint64_t>(ifile));

    const size_t num_faces = chi::ReadBinaryValue<size_t>(ifile);
    for (size_t f = 0; f < num_faces; ++f)
    {
      UnpartitionedMesh::LightWeightFace new_face;
      const size_t num_face_vids = chi::ReadBinaryValue<size_t>(ifile);
      for (size_t v = 0; v < num_face_vids; ++v)
        new_face.vertex_ids.push_back(chi::ReadBinaryValue<uint64_t>(ifile));

      new_face.has_neighbor = chi::ReadBinaryValue<bool>(ifile);
      new_face.neighbor = chi::ReadBinaryValue<uint64_t>(ifile);

      new_cell.faces.push_back(std::move(new_face));
    } // for f

    cells.insert(
      std::make_pair(CellPIDGID(cell_pid, cell_gid), std::move(new_cell)));
  } // for cell c

  //================================================== Read the vertices
  for (size_t v = 0; v < num_vertices; ++v)
  {
    const uint64_t vid = chi::ReadBinaryValue<uint64_t>(ifile);
    const chi_mesh::Vector3 vertex =
      chi::ReadBinaryValue<chi_mesh::Vector3>(ifile);
    vertices.insert(std::make_pair(vid, vertex));
  } // for vertex v

  ifile.close();

  return info_block;
}

std::shared_ptr<MeshContinuum>
SplitFileMeshGenerator::SetupLocalMesh(SplitMeshInfo& mesh_info)
{
  auto grid_ptr = chi_mesh::MeshContinuum::New();

  grid_ptr->GetBoundaryIDMap() = mesh_info.boundary_id_map_;

  auto& cells = mesh_info.cells_;
  auto& vertices = mesh_info.vertices_;

  for (const auto& [vid, vertex] : vertices)
    grid_ptr->vertices.Insert(vid, vertex);

  for (const auto& [pidgid, raw_cell] : cells)
  {
    const auto& [cell_pid, cell_global_id] = pidgid;
    auto cell = SetupCell(
      raw_cell, cell_global_id, cell_pid, STLVertexListHelper(vertices));

    grid_ptr->cells.push_back(std::move(cell));
  }

  SetGridAttributes(
    *grid_ptr,
    static_cast<MeshAttributes>(mesh_info.mesh_attributes_),
    {mesh_info.ortho_Nx_, mesh_info.ortho_Ny_, mesh_info.ortho_Nz_});

  grid_ptr->SetGlobalVertexCount(mesh_info.num_global_vertices_);

  ComputeAndPrintStats(*grid_ptr);

  return grid_ptr;
}

} // namespace chi_mesh