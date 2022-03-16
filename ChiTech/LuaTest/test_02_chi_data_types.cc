#include "unit_tests.h"

#include "ChiDataTypes/byte_array.h"
#include "ChiMesh/Cell/cell.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include "chi_mpi_utils_map_all2all.h"

#include <map>

bool chi_unit_tests::Test_chi_data_types(bool verbose)
{
  bool passed = true;
  std::stringstream output;

  //======================================================= Byte array write/read
  output << "Testing chi_data_types::ByteArray Write and Read\n";
  chi_data_types::ByteArray barr;

  barr.Write<double>(1.01234567890123456789);
  barr.Write<int>(-15600700);
  barr.Write<double>(2.0987654321);
  barr.Write<bool>(false);
  barr.Write<bool>(true);
  barr.Write<chi_mesh::Vector3>(chi_mesh::Vector3(3.0,2.0,1.0));

  auto dbl_value = barr.Read<double>();
  auto int_value = barr.Read<int>();
  auto dbl_value2= barr.Read<double>();
  auto bl_value1 = barr.Read<bool>();
  auto bl_value2 = barr.Read<bool>();
  auto vec3       = barr.Read<chi_mesh::Vector3>();

  if (dbl_value != 1.01234567890123456789 or
      int_value != -15600700 or
      dbl_value2 != 2.0987654321 or
      bl_value1 or
      !bl_value2 or
      vec3[0] != chi_mesh::Vector3(3.0,2.0,1.0)[0] or
      vec3[1] != chi_mesh::Vector3(3.0,2.0,1.0)[1] or
      vec3[2] != chi_mesh::Vector3(3.0,2.0,1.0)[2])

  {
    passed = false;
    output << std::string("chi_data_types::ByteArray Write/Read ... Failed\n");
  }
  else
    output << std::string("chi_data_types::ByteArray Write/Read ... Passed\n");

  if (verbose)
    chi_log.Log() << output.str();

  output.str(""); //clear the stream

  //======================================================= Testing Byte array serialization
  output << "Testing chi_data_types::ByteArray Serialization/DeSerialization\n";
  if (chi_mpi.process_count < 2)
    throw std::logic_error("chi_unit_tests::Test_chi_data_types requires at least"
                           " 2 MPI processes.");

  std::map<int/*pid*/, chi_data_types::ByteArray> send_data;

  chi_mesh::Cell poster_child_cell(chi_mesh::CellType::POLYHEDRON,
                                   chi_mesh::CellType::HEXAHEDRON);
  {
    poster_child_cell.global_id    = 321;
    poster_child_cell.local_id     = 123;
    poster_child_cell.partition_id = 0;
    poster_child_cell.centroid     = chi_mesh::Vector3(0.5, 0.5, 0.5);
    poster_child_cell.material_id  = -2;

    poster_child_cell.vertex_ids = {0, 1, 2, 3,
                                    4, 5, 6, 7};

    //Bottom face
    {
      chi_mesh::CellFace face;
      face.vertex_ids   = {0,3,2,1};
      face.normal       = {0,0,-1};
      face.centroid     = {0.5,0.5,0.0};
      face.has_neighbor = false;
      face.neighbor_id  = 0;
      poster_child_cell.faces.push_back(std::move(face));
    }
    //Top face
    {
      chi_mesh::CellFace face;
      face.vertex_ids   = {4,5,6,7};
      face.normal       = {0,0,1};
      face.centroid     = {0.5,0.5,1.0};
      face.has_neighbor = false;
      face.neighbor_id  = 1;
      poster_child_cell.faces.push_back(std::move(face));
    }
    //Left face
    {
      chi_mesh::CellFace face;
      face.vertex_ids   = {0,4,7,3};
      face.normal       = {-1,0,0};
      face.centroid     = {0.0,0.5,0.5};
      face.has_neighbor = false;
      face.neighbor_id  = 2;
      poster_child_cell.faces.push_back(std::move(face));
    }
    //Right face
    {
      chi_mesh::CellFace face;
      face.vertex_ids   = {1,2,6,5};
      face.normal       = {1,0,0};
      face.centroid     = {1.0,0.5,0.5};
      face.has_neighbor = false;
      face.neighbor_id  = 3;
      poster_child_cell.faces.push_back(std::move(face));
    }
    //Front face
    {
      chi_mesh::CellFace face;
      face.vertex_ids   = {0,1,5,4};
      face.normal       = {0,-1,0};
      face.centroid     = {0.5,0.0,0.5};
      face.has_neighbor = false;
      face.neighbor_id  = 4;
      poster_child_cell.faces.push_back(std::move(face));
    }
    //Back face
    {
      chi_mesh::CellFace face;
      face.vertex_ids   = {3,7,6,2};
      face.normal       = {0,1,0};
      face.centroid     = {0.5,1.0,0.5};
      face.has_neighbor = false;
      face.neighbor_id  = 5;
      poster_child_cell.faces.push_back(std::move(face));
    }
  }

  if (chi_mpi.location_id == 0)
  {
    send_data[1].Append(poster_child_cell.Serialize());
    send_data[1].Append(poster_child_cell.Serialize());
  }

  std::map<int/*pid*/, std::vector<std::byte>> send_data_bytes;

  for (const auto& pid_byte_array : send_data)
    send_data_bytes[pid_byte_array.first] = pid_byte_array.second.Data();

  std::map<int/*pid*/, std::vector<std::byte>> recv_data_bytes =
    chi_mpi_utils::MapAllToAll(send_data_bytes, MPI_BYTE);

  for (const auto& pid_vec_bytes : recv_data_bytes)
  {
    auto& pid = pid_vec_bytes.first;
    auto& vec_bytes = pid_vec_bytes.second;

    chi_data_types::ByteArray byte_array(vec_bytes);

    size_t address = 0;
    while (address < byte_array.Size())
    {
      const chi_mesh::Cell read_cell = chi_mesh::Cell::DeSerialize(byte_array,address);

      auto& rcell = read_cell;
      auto& pcell = poster_child_cell;

      if (rcell.Type()       != pcell.Type()      ) {passed = false; output << "Line: " << __LINE__ << "\n"; break;}
      if (rcell.SubType()    != pcell.SubType()   ) {passed = false; output << "Line: " << __LINE__ << "\n"; break;}
      if (rcell.global_id    != pcell.global_id   ) {passed = false; output << "Line: " << __LINE__ << "\n"; break;}
      if (rcell.local_id     != pcell.local_id    ) {passed = false; output << "Line: " << __LINE__ << "\n"; break;}
      if (rcell.partition_id != pcell.partition_id) {passed = false; output << "Line: " << __LINE__ << "\n"; break;}
      if (rcell.material_id  != pcell.material_id ) {passed = false; output << "Line: " << __LINE__ << "\n"; break;}
      if (rcell.vertex_ids   != pcell.vertex_ids  ) {passed = false; output << "Line: " << __LINE__ << "\n"; break;}

      if (rcell.faces.size() != pcell.faces.size()) {passed = false; output << "Line: " << __LINE__ << "\n"; break;}

      size_t f = 0;
      for (const auto& rface : rcell.faces)
      {
        const auto& pface = pcell.faces[f];

        if (rface.vertex_ids   != pface.vertex_ids)   {passed = false; output << "Line: " << __LINE__ << "\n"; break;}
        if (rface.has_neighbor != pface.has_neighbor) {passed = false; output << "Line: " << __LINE__ << "\n"; break;}
        if (rface.neighbor_id  != pface.neighbor_id)  {passed = false; output << "Line: " << __LINE__ << "\n"; break;}
        ++f;
      }
    }
  }

  if (not passed)
  {
    output << "chi_data_types::ByteArray Serialization/DeSerialization ... Failed\n";
  }
  else
    output << "chi_data_types::ByteArray Serialization/DeSerialization ... Passed\n";


  if (verbose)
    chi_log.Log(LOG_ALL) << output.str();

  return passed;
}