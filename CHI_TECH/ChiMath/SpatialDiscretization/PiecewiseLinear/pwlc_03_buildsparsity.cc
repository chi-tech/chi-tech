#include "pwlc.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include <algorithm>

//###################################################################
/**Builds the sparsity pattern for a Continuous Finite Element Method.*/
void SpatialDiscretization_PWLC::
BuildSparsityPattern(chi_mesh::MeshContinuumPtr grid,
                     std::vector<int> &nodal_nnz_in_diag,
                     std::vector<int> &nodal_nnz_off_diag,
                     chi_math::UnknownManager* unknown_manager)
{
  //======================================== Determine global domain ownership
  std::vector<int> locI_block_addr(chi_mpi.process_count, 0);
  MPI_Allgather(&local_block_address, 1, MPI_INT,
                locI_block_addr.data()   , 1, MPI_INT,
                MPI_COMM_WORLD);

  if (chi_mpi.location_id == 0)
    for (auto locI : locI_block_addr)
      chi_log.Log(LOG_ALLVERBOSE_1) << "Block address = " << locI;
  MPI_Barrier(MPI_COMM_WORLD);

  //**************************************** DEFINE UTILITIES

  //=================================== Dof-handler
  struct DOFHandler
  {
    int local_block_start = 0;
    int local_block_end = 0;
    std::vector<int> locI_block_addr;

    bool IsMapLocal(int ir)
    {return ( ir>=local_block_start and ir<local_block_end);}

    int MapIRLocal(int ir)
    {return ir - local_block_start;}

    int GetLocFromIR(int ir)
    {
      int locI = std::upper_bound(locI_block_addr.begin(),
                                  locI_block_addr.end(),
                                  ir) - locI_block_addr.begin() - 1;
      return locI;
    }

  } dof_handler;

  //=================================== Set dof handler values
  dof_handler.local_block_start = local_block_address;
  dof_handler.local_block_end   = local_block_address +
                                  local_base_block_size;
  dof_handler.locI_block_addr = locI_block_addr;

  // Writes a message on ir error
  auto IR_MAP_ERROR = [] ()
  {
    chi_log.Log(LOG_ALLERROR)
      << "PWL-MapCFEMDOF: ir Mapping error node ";
    exit(EXIT_FAILURE);
  };

  // Writes a message on jr error
  auto JR_MAP_ERROR = [] ()
  {
    chi_log.Log(LOG_ALLERROR)
      << "PWL-MapCFEMDOF: jr Mapping error node ";
    exit(EXIT_FAILURE);
  };

  // Checks whether an integer is already in a vector
  auto IS_VALUE_IN_VECTOR = [](const std::vector<int>& vec, int val)
  {
    bool already_there = false;
    for (auto check_val : vec)
      if (check_val == val) {already_there = true; break;}
    return already_there;
  };

  //**************************************** END OF UTILITIES



  //======================================== Build local sparsity pattern
  chi_log.Log(LOG_0VERBOSE_1) << "Building local sparsity pattern.";
  int local_node_count = local_base_block_size;
  std::vector<std::vector<int>> nodal_connections(local_node_count);

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag .resize(local_node_count, 0);
  nodal_nnz_off_diag.resize(local_node_count, 0);

  for (auto& cell : grid->local_cells)
  {
    for (auto ivid : cell.vertex_ids)
    {
      int ir = MapDOF(ivid); if (ir < 0) IR_MAP_ERROR();

      if (dof_handler.IsMapLocal(ir))
      {
        int il = dof_handler.MapIRLocal(ir);
        std::vector<int>& node_links = nodal_connections[il];

        for (auto& jvid : cell.vertex_ids)
        {
          int jr = MapDOF(jvid); if (jr < 0) JR_MAP_ERROR();

          if (IS_VALUE_IN_VECTOR(node_links,jr)) continue;

          node_links.push_back(jr);
          if (dof_handler.IsMapLocal(jr))
            nodal_nnz_in_diag[il]+=1;
          else
            nodal_nnz_off_diag[il]+=1;
        }//for j
      }//if i local
    }//for i
  }//for cell




  //======================================== Build non-local sparsity pattern
  chi_log.Log(LOG_0VERBOSE_1) << "Building non-local sparsity pattern.";

  // In this process we build a list
  // of ir-nodes that are not local. Each ir-node needs to
  // be furnished with the jr-nodes it links to.

  typedef std::pair<int,std::vector<int>> ROWJLINKS;
  std::vector<ROWJLINKS> ir_links;

  for (auto& cell : grid->local_cells)
  {
    for (auto ivid : cell.vertex_ids)
    {
      int ir = MapDOF(ivid); if (ir < 0) IR_MAP_ERROR();

      if (not dof_handler.IsMapLocal(ir))
      {
        ROWJLINKS  new_ir_link;
        ROWJLINKS* cur_ir_link = &new_ir_link;

        //============================= Check if ir already there
        bool ir_already_there = false;
        for (auto& ir_link : ir_links)
          if (ir == ir_link.first){
            ir_already_there = true;
            cur_ir_link = &ir_link;
            break;
          }

        //============================= Now add links
        auto& node_links = cur_ir_link->second;
        for (auto& jvid : cell.vertex_ids)
        {
          int jr = MapDOF(jvid); if (jr < 0) JR_MAP_ERROR();

          if (IS_VALUE_IN_VECTOR(node_links,jr)) continue;
          else
            node_links.push_back(jr);
        }//for j

        //============================= If its not add it
        if (not ir_already_there)
        {
          cur_ir_link->first = ir;
          ir_links.push_back(*cur_ir_link);
        }

      }//if i not local
    }//for i
  }//for cell

  //======================================== Build communication structure
  chi_log.Log(LOG_0VERBOSE_1) << "Building communication structure.";

  //=================================== Step 1
  // We now serialize the non-local data
  std::vector<std::vector<int>> locI_serialized(chi_mpi.process_count);

  for (const auto& ir_linkage : ir_links)
  {
    int locI = dof_handler.GetLocFromIR(ir_linkage.first);

    locI_serialized[locI].push_back(ir_linkage.second.size()); //row cols amount
    locI_serialized[locI].push_back(ir_linkage.first);         //row num
    for (int jr : ir_linkage.second)
      locI_serialized[locI].push_back(jr);                     //col num
  }

  //=================================== Step 2
  // Establish the size of the serialized data
  // to send to each location and communicate
  // to get receive count.
  std::vector<int> sendcount(chi_mpi.process_count, 0);
  std::vector<int> recvcount(chi_mpi.process_count, 0);
  int locI=0;
  for (const auto& locI_data : locI_serialized)
  {
    sendcount[locI] = locI_data.size();

    if (chi_mpi.location_id == 0)
      chi_log.Log(LOG_ALLVERBOSE_1)
        << "To send to " << locI
        << " = " << sendcount[locI];

    ++locI;
  }

  MPI_Alltoall(sendcount.data(), 1, MPI_INT,
               recvcount.data(), 1, MPI_INT,
               MPI_COMM_WORLD);

  //=================================== Step 3
  // We now establish send displacements and
  // receive displacements.
  std::vector<int> send_displs(chi_mpi.process_count,0);
  std::vector<int> recv_displs(chi_mpi.process_count,0);

  int send_displ_c = 0;
  int recv_displ_c = 0;

  int c=0;
  for (int send_count : sendcount)
  {
    send_displs[c++] = send_displ_c;
    send_displ_c += send_count;
  }
  c=0;
  for (int recv_count : recvcount)
  {
    recv_displs[c++] = recv_displ_c;
    recv_displ_c += recv_count;
  }

  //======================================== Communicate data
  chi_log.Log(LOG_0VERBOSE_1) << "Communicating non-local rows.";

  // We now initialize the buffers and
  // communicate the data
  std::vector<int> sendbuf;
  std::vector<int> recvbuf;

  sendbuf.reserve(send_displ_c);
  recvbuf.resize(recv_displ_c,0);

  for (const auto& serial_block : locI_serialized)
    for (int data_val : serial_block)
      sendbuf.push_back(data_val);

  MPI_Alltoallv(sendbuf.data(),
                sendcount.data(),
                send_displs.data(),
                MPI_INT,
                recvbuf.data(),
                recvcount.data(),
                recv_displs.data(),
                MPI_INT,
                MPI_COMM_WORLD);

  //======================================== Deserialze data
  chi_log.Log(LOG_0VERBOSE_1) << "Deserialize data.";

  std::vector<ROWJLINKS> foreign_ir_links;

  for (size_t k=0; k<recvbuf.size(); )
  {
    int num_values = recvbuf[k++];
    int ir         = recvbuf[k++];

    ROWJLINKS new_links;
    new_links.first = ir;
    new_links.second.reserve(num_values);
    for (int i=0; i<num_values; ++i)
      new_links.second.push_back(recvbuf[k++]);

    foreign_ir_links.push_back(std::move(new_links));
  }

  //======================================== Adding to sparsity pattern
  for (const auto& ir_linkage : foreign_ir_links)
  {
    int ir = ir_linkage.first;

    if (not dof_handler.IsMapLocal(ir)) IR_MAP_ERROR();

    int il = dof_handler.MapIRLocal(ir);
    std::vector<int>& node_links = nodal_connections[il];

    for (int jr : ir_linkage.second)
    {
      if (IS_VALUE_IN_VECTOR(node_links,jr)) continue;

      node_links.push_back(jr);
      if (dof_handler.IsMapLocal(jr))
        nodal_nnz_in_diag[il]+=1;
      else
        nodal_nnz_off_diag[il]+=1;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //======================================== Spacing according to unknown
  //                                         manager
  if (unknown_manager != nullptr)
  {
    auto backup_nnz_in_diag  = nodal_nnz_in_diag;
    auto backup_nnz_off_diag = nodal_nnz_off_diag;

    unsigned int N = unknown_manager->GetTotalUnknownStructureSize();

    nodal_nnz_in_diag.clear();
    nodal_nnz_off_diag.clear();

    nodal_nnz_in_diag.resize(local_base_block_size*N,0);
    nodal_nnz_off_diag.resize(local_base_block_size*N,0);

    if (unknown_manager->dof_storage_type == chi_math::UnknownStorageType::NODAL)
    {
      int ir = -1;
      for (int i=0; i<local_base_block_size; ++i)
      {
        for (int j=0; j<N; ++j)
        {
          ++ir;
          nodal_nnz_in_diag[ir] = backup_nnz_in_diag[i];
          nodal_nnz_off_diag[ir] = backup_nnz_off_diag[i];
        }//for j
      }//for i
    }
    else if (unknown_manager->dof_storage_type == chi_math::UnknownStorageType::BLOCK)
    {
      int ir = -1;
      for (int j=0; j<N; ++j)
      {
        for (int i=0; i<local_base_block_size; ++i)
        {
          ++ir;
          nodal_nnz_in_diag[ir] = backup_nnz_in_diag[i];
          nodal_nnz_off_diag[ir] = backup_nnz_off_diag[i];
        }//for i
      }//for j
    }

  }//if unknown manager supplied
}