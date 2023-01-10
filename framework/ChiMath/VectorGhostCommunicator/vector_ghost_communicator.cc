#include "vector_ghost_communicator.h"

#include <map>
#include <string>

#include "ChiMPI/chi_mpi_utils_map_all2all.h"

chi_math::VectorGhostCommunicator::
  VectorGhostCommunicator(uint64_t local_size,
                          uint64_t global_size,
                          std::vector<int64_t> ghost_indices,
                          MPI_Comm communicator) :
  m_local_size(local_size),
  m_globl_size(global_size),
  m_ghost_indices(std::move(ghost_indices)),
  m_comm(communicator)
{
  const std::string fname = "chi_math::VectorWithGhosts::VectorWithGhosts";
  MPI_Comm_rank(m_comm, &m_location_id);
  MPI_Comm_size(m_comm, &m_process_count);

  std::vector<uint64_t> locI_local_size(m_process_count, 0);

  MPI_Allgather(&m_local_size,          //sendbuf
                1, MPI_UINT64_T,        //sendcount + sendtype
                locI_local_size.data(), //recvbuf
                1, MPI_UINT64_T,        //recvcount + recvtype
                m_comm);                //communicator

  //======================================== Determine location extents
  //This will allow us to locally determine
  //entry ownership
  m_locI_extents.assign(m_process_count+1,0);
  uint64_t offset = 0;
  for (size_t locI=0; locI<m_process_count; ++locI)
  {
    m_locI_extents[locI] = offset;
    offset += locI_local_size[locI];
  }
  m_locI_extents.back() = offset;

  //======================================== READ THIS
  // The next steps can be a little confusing. What we essentially want to
  // know is the information necessary to perform a vector-scatter operation
  // multiple times and as efficiently as possible. This information includes
  // what information to send to other locations AS WELL AS what information
  // we will be receiving.
  // Definitions:
  //   giver - locations that need to send this location something
  //   taker - locations to whom this location needs to send something

  //======================================== Build giver query-map
  //This map is keyed by giver location.
  std::map<int, std::vector<int64_t>> giver_query_map;
  for (int64_t ggid : m_ghost_indices)
    giver_query_map[FindOwnerPID(ggid)].push_back(ggid);

  //======================================== Build recv-map
  //Buy building the query-map we essentially
  //the order in which the ghost values will be
  //received. Therefore, we have to keep track
  //of these changes by storing a map of their
  //reconfiguration.
  //This map will be used during a communication.
  size_t recv_id=0;
  for (const auto& pid_list : giver_query_map)
    for (int64_t gid : pid_list.second)
    {
      m_ghost_ids_to_recv_map[gid] = recv_id;
      ++recv_id;
    }

  //======================================== Build recvcounts and recvdispls
  //Since the giver query-map contains the ids of information to be obtained
  //from the givers we can also build recvcounts and recvdispls from it
  //since, during each scatter, this is the configuration in which the info
  //will come from the givers. This will be used in a call to
  //MPI_Alltoallv.
  m_recvcounts.assign(m_process_count, 0);
  m_recvdispls.assign(m_process_count, 0);
  uint64_t total_recvcounts = 0;
  for (const auto& pid_vec_pair : giver_query_map)
  {
    const int    pid        = pid_vec_pair.first;
    const size_t num_counts = pid_vec_pair.second.size();
    m_recvcounts[pid] = static_cast<int>(num_counts);
    m_recvdispls[pid] = static_cast<int>(total_recvcounts);
    total_recvcounts += num_counts;
  }

  //======================================== Communicate giver query map
  //Here the giver_query_map is transformed into
  //a taker_query_map. This is now a different perspective
  //and this location now has a map of what data is needed
  //by each dependent process.
  auto taker_query_map = chi_mpi_utils::MapAllToAll(
    giver_query_map, //map key+list
    MPI_INT64_T,     //datatype of list
    m_comm           //communicator
  );

  //======================================== Map the query indices to
  //                                         local indices
  //First determine how much data we will
  //eventually have to send. This allows us
  //to reserve some space.
  size_t amount_to_give = 0;
  for (const auto& pid_list_pair : taker_query_map)
    amount_to_give += pid_list_pair.second.size();

  //======================================== Premap local ids
  //These will be used during a communication.
  m_takers_local_ids.clear();
  m_takers_local_ids.reserve(amount_to_give);
  for (const auto& pid_list_pair : taker_query_map)
  {
    const auto& list = pid_list_pair.second;
    for (int64_t gid : list)
    {
      const uint64_t local_block_beg = m_locI_extents[m_location_id];
      const uint64_t local_block_end = m_locI_extents[m_location_id+1];
      if (not (gid >= local_block_beg and gid <  local_block_end))
        throw std::logic_error(fname + ": Error mapping index in "
                                       "takers_query_map.");

      m_takers_local_ids.push_back(gid - static_cast<int64_t>(local_block_beg));
    }//for gid
  }

  //======================================== Build sendcounts and senddispls
  //Since we now know what the structure will
  //of sending the data we can compute m_sendcounts
  //and m_senddispls, which will be used in a call
  //to MPI_Alltoallv
  m_sendcounts.assign(m_process_count, 0);
  m_senddispls.assign(m_process_count, 0);

  uint64_t total_sendcounts = 0;
  for (const auto& pid_vec_pair : taker_query_map)
  {
    const int    pid        = pid_vec_pair.first;
    const size_t num_counts = pid_vec_pair.second.size();
    m_sendcounts[pid] = static_cast<int>(num_counts);
    m_senddispls[pid] = static_cast<int>(total_sendcounts);
    total_sendcounts += num_counts;
  }
}


//###################################################################
int chi_math::VectorGhostCommunicator::FindOwnerPID(uint64_t global_id) const
{
  if (global_id >= m_globl_size)
    throw std::logic_error(
      "chi_math::VectorWithGhosts:FindOwnerPID global_id >= m_globl_size");
  for (int locI=0; locI<static_cast<int>(m_process_count); ++locI)
    if (global_id >= m_locI_extents[locI] and
        global_id < m_locI_extents[locI+1])
      return locI;
  return int(-1);
}

//###################################################################
void chi_math::VectorGhostCommunicator::
  CommunicateGhostEntries(std::vector<double> &local_vector) const
{
  if (local_vector.size() != (m_local_size + m_ghost_indices.size()))
    throw std::logic_error(
      "chi_math::VectorWithGhosts::CommunicateGhostEntries: Vector size "
      "mismatch.");

  //======================================== Build data to be sent
  const size_t amount_to_send = m_takers_local_ids.size();
  std::vector<double> send_data;
  send_data.reserve(amount_to_send);
  for (const int64_t local_id : m_takers_local_ids)
    send_data.push_back(local_vector[local_id]);

  //======================================== Allocate receive buffer
  const size_t amount_to_recv = m_ghost_indices.size();
  std::vector<double> recv_data(amount_to_recv, 0.0);

  //======================================== Communicate
  MPI_Alltoallv(send_data.data(),
                m_sendcounts.data(),
                m_senddispls.data(),
                MPI_DOUBLE,
                recv_data.data(),
                m_recvcounts.data(),
                m_recvdispls.data(),
                MPI_DOUBLE,
                m_comm);

  //======================================== Populate local vector with ghost
  //                                         data
  for (size_t k=0; k<amount_to_recv; ++k)
  {
    const int64_t gid = m_ghost_indices[k];
    const size_t  rid = m_ghost_ids_to_recv_map.at(gid);
    local_vector[m_local_size + k] = recv_data[rid];
  }
}


//#########################################################
/**Returns a vector with enough space to hold both `local_size` entries,
 * as specified in the VectorGhostCommunicator-constructor, and entries
 * associated with the ghost-ids. The entries default to zero.*/
std::vector<double> chi_math::VectorGhostCommunicator::
  MakeGhostedVector() const
{
  std::vector<double> vec(m_local_size+m_ghost_indices.size(),0.0);
  return vec;
}


//#########################################################
/**Returns a vector with enough space to hold both `local_size` entries,
 * as specified in the VectorGhostCommunicator-constructor, and entries
 * associated with the ghost-ids. The local entries are copied from the
 * vector supplied as an argument whilst the ghosts are defaulted to zero.*/
std::vector<double> chi_math::VectorGhostCommunicator::
  MakeGhostedVector(const std::vector<double>& unghosted_vector) const
{
  if (unghosted_vector.size() < m_local_size)
    throw std::logic_error(
      "chi_math::VectorGhostCommunicator::MakeGhostedVector: "
      "unghosted_vector.size() < m_local_size");

  std::vector<double> vec(m_local_size+m_ghost_indices.size(),0.0);
  for (size_t i=0; i<m_local_size; ++i)
    vec[i] = unghosted_vector[i];
  return vec;
}