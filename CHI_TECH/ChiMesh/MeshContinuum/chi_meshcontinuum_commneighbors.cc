#include "chi_meshcontinuum.h"

#include "ChiMesh/Cell/cell_slab.h"
#include "ChiMesh/Cell/cell_polygon.h"
#include "ChiMesh/Cell/cell_polyhedron.h"

#include <chi_log.h>
#include <chi_mpi.h>
extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Communicates neighboring cells to this location for use by methods
 * such as the interior penalty method. The method populates the
 * supplied vector neighbor_cells. The complete cell is
 * not populated, the face neighbors are not transferred.*/
void chi_mesh::MeshContinuum::CommunicatePartitionNeighborCells(
  std::vector<chi_mesh::Cell*>& neighbor_cells)
{
  chi_log.Log(LOG_0)
    << "Communicating partition neighbors.";
  MPI_Barrier(MPI_COMM_WORLD);

  std::set<int> local_neighboring_cell_indices;
  std::set<int> neighboring_partitions;

  //============================================= Collect local neighboring
  //                                              cell indices and neighboring
  //                                              partition indices
  // First establish which local nodes are
  // neighbors to other partitions. Build a set
  // of local cells and a set of destination
  // partitions.
  for (auto& cell : local_cells)
  {
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor and (not face.IsNeighborLocal(*this)))
      {
        local_neighboring_cell_indices.insert(cell.local_id);
        neighboring_partitions.insert(cells[face.neighbor_id].partition_id);
      }
    }
  }

  //============================================= Subscribe neighbor-partitions
  //                                              the local cells they need
  // Now we establish a pair for each unique location.
  // pair.first is the partition-id of the location.
  // pair.second is the list of cells to be sent to that location.
  size_t num_neighbor_partitions = neighboring_partitions.size();

  typedef std::pair<int,std::vector<int>> ListOfCells;

  std::vector<ListOfCells> destination_subscriptions;
  destination_subscriptions.reserve(num_neighbor_partitions);
  for (int adj_part : neighboring_partitions)
  {
    ListOfCells new_list;

    new_list.first = adj_part;
    for (int local_cell_index : local_neighboring_cell_indices)
    {
      auto& cell = local_cells[local_cell_index];

      for (auto& face : cell.faces)
      {
        if ((face.has_neighbor) and (not face.IsNeighborLocal(*this)) )
        {
          if (cells[face.neighbor_id].partition_id == adj_part)
            new_list.second.push_back(local_cell_index);
        }
      }//for faces
    }//for neighbor cells

    destination_subscriptions.push_back(new_list);
  }//for adj partition

  //============================================= Serialize
  // For each location we now serialize the cells
  // it needs.
  // The serialized values will be as follows
  // - cell_type
  // - cell_glob_index

  // - cell_local_index

  // - cell_mat_id
  // - cell_dof_count
  // - cell_face_count
  //
  // - dof 0 glob_index
  //     to
  // - dof N glob_index
  //
  // - face_0 dof_count
  // - face_0 dof 0 glob_index
  //     to
  // - face_0 dof fN glob_index
  //
  // - repeat all face info
  std::vector<ListOfCells> destination_serialized_data;

  destination_serialized_data.reserve(num_neighbor_partitions);

  int total_serial_size = 0;
  for (auto& cell_list : destination_subscriptions)
  {
    ListOfCells new_serial_data;

    new_serial_data.first = cell_list.first;
    for (int local_cell_index : cell_list.second)
    {
      auto& cell = local_cells[local_cell_index];

      std::vector<int>& border_cell_info = new_serial_data.second;

      if (cell.Type() == chi_mesh::CellType::SLAB)
        border_cell_info.push_back(3);                         //cell_type
      else if (cell.Type() == chi_mesh::CellType::POLYGON)
        border_cell_info.push_back(4);                         //cell_type
      else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
        border_cell_info.push_back(5);                         //cell_type
      else
        border_cell_info.push_back(-1);                        //cell_type

      border_cell_info.push_back(cell.global_id);         //cell_glob_index
      border_cell_info.push_back(cell.local_id);         //cell_local_index
      border_cell_info.push_back(cell.material_id);            //cell_mat_id
      border_cell_info.push_back(cell.partition_id);           //cell_par_id
      border_cell_info.push_back(cell.vertex_ids.size());      //cell_dof_count
      border_cell_info.push_back(cell.faces.size());           //cell_face_count

      for (auto vid : cell.vertex_ids) //vid = vertex-id
        border_cell_info.push_back(vid);//dof 0 to N

      for (auto& face : cell.faces)
      {
        int face_dof_count = face.vertex_ids.size();
        border_cell_info.push_back(face_dof_count);         //face dof_count
        for (int fvid : face.vertex_ids) //fvid = face-vertex-id
          border_cell_info.push_back(fvid);
        //face dof 0 to fN
      }
    }
    total_serial_size += new_serial_data.second.size();
    destination_serialized_data.push_back(new_serial_data);
  }

  //============================================= Build send-counts and
  //                                              send-displacements arrays, and
  //                                              serialized vector
  std::vector<int> send_counts(chi_mpi.process_count,0);
  std::vector<int> send_displs(chi_mpi.process_count,0);
  std::vector<int> global_serialized_data;

  global_serialized_data.reserve(total_serial_size);

  int displacement = 0;
  for (auto& info : destination_serialized_data)
  {
    send_counts[info.first] = info.second.size();
    send_displs[info.first] = displacement;
    displacement += info.second.size();

    for (int val : info.second)
      global_serialized_data.push_back(val);
  }

  //============================================= Communicate counts
  std::vector<int> recv_counts(chi_mpi.process_count,0);

  MPI_Alltoall(send_counts.data(), 1, MPI_INT,
               recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  //============================================= Build receive displacements
  std::vector<int> recv_displs(chi_mpi.process_count,0);
  int total_receive_size = 0;
  int c=0;
  for (auto val : recv_counts)
  {
    recv_displs[c] = total_receive_size;
    total_receive_size += val;
    ++c;
  }

  //============================================= Receive serialized data
  std::vector<int> global_receive_data(total_receive_size,0);
  MPI_Alltoallv(global_serialized_data.data(),
                send_counts.data(),
                send_displs.data(),
                MPI_INT,
                global_receive_data.data(),
                recv_counts.data(),
                recv_displs.data(),
                MPI_INT,
                MPI_COMM_WORLD);

  //============================================= Deserialize
  {
    int k=0;
    while (k<global_receive_data.size())
    {

      int cell_type       = global_receive_data[k]; k++;
      chi_mesh::Cell* cell;
      if (cell_type == 3)
        cell = new chi_mesh::CellSlab;
      else if (cell_type == 4)
        cell = new chi_mesh::CellPolygon;
      else if (cell_type == 5)
        cell = new chi_mesh::CellPolyhedron;
      else
      {
        chi_log.Log(LOG_ALLERROR)
          << "chi_mesh::MeshContinuum::CommunicatePartitionNeighborCells, "
          << "unsupported cell type encountered during deserialization."
          << " " << cell_type << " " << k << " " << global_receive_data.size();
        exit(EXIT_FAILURE);
      }

      cell->global_id = global_receive_data[k]; k++;
      cell->local_id = global_receive_data[k]; k++;
      cell->material_id    = global_receive_data[k]; k++;
      cell->partition_id   = global_receive_data[k]; k++;
      int cell_dof_count  = global_receive_data[k]; k++;
      int cell_face_count = global_receive_data[k]; k++;

      cell->vertex_ids.reserve(cell_dof_count);
      for (int v=0; v<cell_dof_count; ++v)
      {
        cell->vertex_ids.push_back(global_receive_data[k]);
        k++;
      }

      cell->centroid = ComputeCentroidFromListOfNodes(cell->vertex_ids);

      for (int f=0; f<cell_face_count; ++f)
      {
        int face_dof_count = global_receive_data[k]; k++;
        cell->faces.emplace_back();
        cell->faces[f].vertex_ids.reserve(face_dof_count);
        for (int fv=0; fv<face_dof_count; fv++)
        {
          int vgi = global_receive_data[k]; k++;
          cell->faces[f].vertex_ids.push_back(vgi);
        }
        cell->faces[f].centroid =
          ComputeCentroidFromListOfNodes(cell->faces[f].vertex_ids);

        if (cell->Type() == chi_mesh::CellType::SLAB)
        {
          if (f == 0)
            cell->faces[f].normal = chi_mesh::Normal(0.0,0.0,-1.0);
          else
            cell->faces[f].normal = chi_mesh::Normal(0.0,0.0, 1.0);
        }
        else if (cell->Type() == chi_mesh::CellType::POLYGON)
        {
          auto k_hat = chi_mesh::Normal(0.0,0.0, 1.0);

          auto& v0 = *vertices[cell->faces[f].vertex_ids[0]];
          auto& v1 = *vertices[cell->faces[f].vertex_ids[1]];

          auto v01 = v1 - v0;

          cell->faces[f].normal = (v01.Cross(k_hat)).Normalized();
        }
        else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
        {
          auto& v0 = *vertices[cell->faces[f].vertex_ids[0]];
          auto& v1 = *vertices[cell->faces[f].vertex_ids[1]];
          auto& v2 = cell->faces[f].centroid;

          auto v01 = v1 - v0;
          auto v12 = v2 - v1;

          cell->faces[f].normal = (v01.Cross(v12)).Normalized();
        }
      }

      neighbor_cells.push_back(cell);
    }
  }



}