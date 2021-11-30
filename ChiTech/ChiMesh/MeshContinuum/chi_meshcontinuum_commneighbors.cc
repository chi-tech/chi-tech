#include "chi_meshcontinuum.h"

#include "ChiMesh/Cell/cell.h"

#include "ChiDataTypes/byte_array.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Communicates neighboring cells to this location for use by methods
 * such as the interior penalty method. The method populates the
 * supplied vector neighbor_cells. The complete cell is
 * not populated, the face neighbors are not transferred.*/
void chi_mesh::MeshContinuum::CommunicatePartitionNeighborCells(
  std::map<uint64_t, chi_mesh::Cell*>& neighbor_cells)
{
  MPI_Barrier(MPI_COMM_WORLD);

  std::set<uint64_t> local_cells_on_partition_bndry;
  std::set<uint64_t> neighboring_partitions;

  //============================================= Collect local neighboring
  //                                              cell indices and neighboring
  //                                              partition indices
  // First establish which local nodes are
  // neighbors to other partitions. Build a set
  // of local cells and a set of destination
  // partitions.
  for (auto& cell : local_cells)
    for (auto& face : cell.faces)
      if (face.has_neighbor and (not face.IsNeighborLocal(*this)))
      {
        local_cells_on_partition_bndry.insert(cell.local_id);
        neighboring_partitions.insert(cells[face.neighbor_id].partition_id);
      }

  //============================================= Subscribe neighbor-partitions
  //                                              the local cells they need
  // Now we establish a pair for each unique location.
  // pair.first is the partition-id of the location.
  // pair.second is the list of cells to be sent to that location.
  size_t num_neighbor_partitions = neighboring_partitions.size();

  typedef std::pair<uint64_t ,std::vector<uint64_t>> LocationListOfCells;

  std::vector<LocationListOfCells> destination_subscriptions;
  destination_subscriptions.reserve(num_neighbor_partitions);
  for (uint64_t adj_part : neighboring_partitions)
  {
    LocationListOfCells new_list;

    new_list.first = adj_part;
    for (uint64_t local_cell_index : local_cells_on_partition_bndry)
    {
      auto& cell = local_cells[local_cell_index];

      for (auto& face : cell.faces)
        if ((face.has_neighbor) and (not face.IsNeighborLocal(*this)) )
          if (cells[face.neighbor_id].partition_id == adj_part)
            new_list.second.push_back(local_cell_index);

    }//for neighbor cells

    destination_subscriptions.push_back(new_list);
  }//for adj partition


  //======================================== Build serial data for specific
  //                                         locations
  // This is needed to later build send-counts
  // and displacements
  std::vector<chi_data_types::ByteArray> locI_serial_data(chi_mpi.process_count);
  for (const auto& cell_list : destination_subscriptions)
  {
    const size_t locID                      = cell_list.first;
    const std::vector<uint64_t>& index_list = cell_list.second;

    for (uint64_t cell_local_id : index_list)
    {
      auto& cell = local_cells[cell_local_id];
      locI_serial_data[locID].Append(cell.Serialize());
    }
  }

  //======================================== Declare communication arrays
  std::vector<int> sendcounts(chi_mpi.process_count, 0);
  std::vector<int> senddispls(chi_mpi.process_count, 0);
  std::vector<int> recvcounts(chi_mpi.process_count, 0);
  std::vector<int> recvdispls(chi_mpi.process_count, 0);

  //======================================== Populate sendcounts and senddispls
  {
    int displacement=0;
    for (int locI=0; locI<chi_mpi.process_count; ++locI)
    {
      sendcounts[locI] = static_cast<int>(locI_serial_data[locI].Size());
      senddispls[locI] = displacement;
      displacement += sendcounts[locI];
    }//for locI
  }

  //======================================== Communicate to populate recvcounts
  MPI_Alltoall(sendcounts.data(), //sendbuf
               1, MPI_INT,        //sendcount, sendtype
               recvcounts.data(), //recvbuf
               1, MPI_INT,        //recvcount, recvtype
               MPI_COMM_WORLD);   //communicator

  //======================================== Populate recvdispls and
  //                                         total_recv_count
  size_t total_recv_count;
  {
    int displacement=0;
    for (int locI=0; locI<chi_mpi.process_count; ++locI)
    {
      recvdispls[locI] = displacement;
      displacement += recvcounts[locI];
    }//for locI
    total_recv_count = displacement;
  }

  //======================================== Make a consolidated send buffer
  chi_data_types::ByteArray sendbuf_serial_data;
  for (const auto& serial_data : locI_serial_data)
    sendbuf_serial_data.Append(serial_data);

  //======================================== Make a recvbuf
  chi_data_types::ByteArray recvbuf(total_recv_count);

  //======================================== Communicate serial data
  MPI_Alltoallv(sendbuf_serial_data.Data().data(), //sendbuf
                sendcounts.data(),                 //sendcounts
                senddispls.data(),                 //senddispls
                MPI_BYTE,                          //sendtype
                recvbuf.Data().data(),             //recvbuf
                recvcounts.data(),                 //recvcounts
                recvdispls.data(),                 //recvdispls
                MPI_BYTE,                          //recvtype
                MPI_COMM_WORLD);                   //comm

  //======================================== Read cells from the received buffer
  if (recvbuf.Size()>0)
  {
    size_t address=0;
    while (address < recvbuf.Size())
    {
      auto cell_read = chi_mesh::Cell::DeSerialize(recvbuf, address);
      auto cell = new chi_mesh::Cell(std::move(cell_read));

      neighbor_cells.insert(std::pair<uint64_t, Cell*>(cell->global_id,cell));
    }
  }

}