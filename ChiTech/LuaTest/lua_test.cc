#include <ChiLua/chi_lua.h>

#include "ChiMath/dynamic_vector.h"
#include "ChiMath/dynamic_matrix.h"

#include "ChiDataTypes/byte_array.h"
#include "ChiMesh/Cell/cell.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**This is a lua test function.
\param argument1 Any Argument of any type.
\ingroup LuaGeneralUtilities
 */
int chiLuaTest(lua_State* L)
{

  chi_math::DynamicVector<double> vec(5, 1.0);
  chi_math::DynamicMatrix<double> mat(5,7,1.0);

  chi_log.Log() << vec.PrintStr() << std::endl;
  chi_log.Log() << "Hello\n" << mat.PrintStr();

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

  chi_log.Log() << dbl_value << " "
                << int_value << " "
                << dbl_value2 << " "
                << bl_value1 << " "
                << bl_value2 << " "
                << vec3.PrintS() << " ";


  std::vector<chi_data_types::ByteArray> locI_local_serialized_cells(chi_mpi.process_count);

  if (chi_mpi.location_id == 0)
  {
    chi_mesh::Cell cell(chi_mesh::CellType::POLYHEDRON,
                        chi_mesh::CellType::HEXAHEDRON);
    cell.global_id    = 321;
    cell.local_id     = 123;
    cell.partition_id = 0;
    cell.centroid     = chi_mesh::Vector3(0.5, 0.5, 0.5);
    cell.material_id  = -2;

    cell.vertex_ids = {0,1,2,3,
                       4,5,6,7};

    //Bottom face
    {
      chi_mesh::CellFace face;
      face.vertex_ids   = {0,3,2,1};
      face.normal       = {0,0,-1};
      face.centroid     = {0.5,0.5,0.0};
      face.has_neighbor = false;
      face.neighbor_id  = 0;
      chi_log.Log() << face.normal.PrintS();
      cell.faces.push_back(std::move(face));
    }
    //Top face
    {
      chi_mesh::CellFace face;
      face.vertex_ids   = {4,5,6,7};
      face.normal       = {0,0,1};
      face.centroid     = {0.5,0.5,1.0};
      face.has_neighbor = false;
      face.neighbor_id  = 1;
      chi_log.Log() << face.normal.PrintS();
      cell.faces.push_back(std::move(face));
    }
    //Left face
    {
      chi_mesh::CellFace face;
      face.vertex_ids   = {0,4,7,3};
      face.normal       = {-1,0,0};
      face.centroid     = {0.0,0.5,0.5};
      face.has_neighbor = false;
      face.neighbor_id  = 2;
      chi_log.Log() << face.normal.PrintS();
      cell.faces.push_back(std::move(face));
    }
    //Right face
    {
      chi_mesh::CellFace face;
      face.vertex_ids   = {1,2,6,5};
      face.normal       = {1,0,0};
      face.centroid     = {1.0,0.5,0.5};
      face.has_neighbor = false;
      face.neighbor_id  = 3;
      chi_log.Log() << face.normal.PrintS();
      cell.faces.push_back(std::move(face));
    }
    //Front face
    {
      chi_mesh::CellFace face;
      face.vertex_ids   = {0,1,5,4};
      face.normal       = {0,-1,0};
      face.centroid     = {0.5,0.0,0.5};
      face.has_neighbor = false;
      face.neighbor_id  = 4;
      chi_log.Log() << face.normal.PrintS();
      cell.faces.push_back(std::move(face));
    }
    //Back face
    {
      chi_mesh::CellFace face;
      face.vertex_ids   = {3,7,6,2};
      face.normal       = {0,1,0};
      face.centroid     = {0.5,1.0,0.5};
      face.has_neighbor = false;
      face.neighbor_id  = 5;
      chi_log.Log() << face.normal.PrintS();
      cell.faces.push_back(std::move(face));
    }

    locI_local_serialized_cells[2].Append(cell.Serialize());
    locI_local_serialized_cells[2].Append(cell.Serialize());
  }

  std::vector<int> locI_sendcounts(chi_mpi.process_count, 0);
  std::vector<int> locI_senddispls(chi_mpi.process_count, 0);
  std::vector<int> locI_recvcounts(chi_mpi.process_count, 0);
  {
    int displacement = 0;
    for (int locI=0; locI<chi_mpi.process_count; ++locI)
    {
      locI_sendcounts[locI] = static_cast<int>(locI_local_serialized_cells[locI].Size());
      locI_senddispls[locI] = displacement;
      displacement += locI_sendcounts[locI];
    }

    MPI_Alltoall(locI_sendcounts.data(), //sendbuf,
                 1, MPI_INT,             //sendcount, sendtype
                 locI_recvcounts.data(), //recvbuf
                 1, MPI_INT,             //recvcount, recvtype
                 MPI_COMM_WORLD);        //comm
  }

  std::vector<int> locI_recvdispls(chi_mpi.process_count, 0);
  size_t total_recv_size=0;
  {
    int displacement = 0;
    for (int locI=0; locI<chi_mpi.process_count; ++locI)
    {
      locI_recvdispls[locI] = displacement;
      displacement += locI_recvcounts[locI];
    }
    total_recv_size = displacement;
  }

  chi_data_types::ByteArray local_aggregated_serial_data;
  std::vector<std::byte>    recvd_aggregated_serial_data(total_recv_size);
  for (int locI=0; locI<chi_mpi.process_count; ++locI)
    local_aggregated_serial_data.Append(locI_local_serialized_cells[locI]);

  MPI_Alltoallv(local_aggregated_serial_data.Data().data(),
                locI_sendcounts.data(),
                locI_senddispls.data(),
                MPI_BYTE,
                recvd_aggregated_serial_data.data(),
                locI_recvcounts.data(),
                locI_recvdispls.data(),
                MPI_BYTE,
                MPI_COMM_WORLD);

  chi_data_types::ByteArray recvd_byte_array;
  recvd_byte_array.Append(recvd_aggregated_serial_data);

  if (recvd_byte_array.Size() > 0)
  {
    size_t address=0;
    while (address < recvd_byte_array.Size())
    {
      chi_log.Log(LOG_ALL) << "Cell at address " << address
                         << " of " << recvd_byte_array.Size();
      chi_mesh::Cell cell = chi_mesh::Cell::DeSerialize(recvd_byte_array,address);
      chi_log.Log(LOG_ALL) << "Address at end " << address;
    }
  }



  return 0;
}