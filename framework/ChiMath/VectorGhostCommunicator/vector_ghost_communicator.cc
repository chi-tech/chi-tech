#include "vector_ghost_communicator.h"

#include <map>
#include <string>

#include "ChiMPI/chi_mpi_utils_map_all2all.h"

chi_math::VectorGhostCommunicator::
  VectorGhostCommunicator(uint64_t local_size,
                          uint64_t global_size,
                          std::vector<int64_t> ghost_indices,
                          MPI_Comm communicator) :
    local_size_(local_size),
    globl_size_(global_size),
    ghost_indices_(std::move(ghost_indices)),
    comm_(communicator)
{
  const std::string fname = "chi_math::VectorWithGhosts::VectorWithGhosts";
  MPI_Comm_rank(comm_, &location_id_);
  MPI_Comm_size(comm_, &process_count_);

  std::vector<uint64_t> locI_local_size(process_count_, 0);

  MPI_Allgather(&local_size_,          //sendbuf
                1, MPI_UINT64_T,        //sendcount + sendtype
                locI_local_size.data(), //recvbuf
                1, MPI_UINT64_T,        //recvcount + recvtype
                comm_);                //communicator

  //======================================== Determine location extents
  //This will allow us to locally determine
  //entry ownership
  locI_extents_.assign(process_count_ + 1, 0);
  uint64_t offset = 0;
  for (size_t locI=0; locI < process_count_; ++locI)
  {
    locI_extents_[locI] = offset;
    offset += locI_local_size[locI];
  }
  locI_extents_.back() = offset;

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
  for (int64_t ggid : ghost_indices_)
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
      ghost_ids_to_recv_map_[gid] = recv_id;
      ++recv_id;
    }

  //======================================== Build recvcounts and recvdispls
  //Since the giver query-map contains the ids of information to be obtained
  //from the givers we can also build recvcounts and recvdispls from it
  //since, during each scatter, this is the configuration in which the info
  //will come from the givers. This will be used in a call to
  //MPI_Alltoallv.
  recvcounts_.assign(process_count_, 0);
  recvdispls_.assign(process_count_, 0);
  uint64_t total_recvcounts = 0;
  for (const auto& pid_vec_pair : giver_query_map)
  {
    const int    pid        = pid_vec_pair.first;
    const size_t num_counts = pid_vec_pair.second.size();
    recvcounts_[pid] = static_cast<int>(num_counts);
    recvdispls_[pid] = static_cast<int>(total_recvcounts);
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
    comm_           //communicator
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
  takers_local_ids_.clear();
  takers_local_ids_.reserve(amount_to_give);
  for (const auto& pid_list_pair : taker_query_map)
  {
    const auto& list = pid_list_pair.second;
    for (int64_t gid : list)
    {
      const uint64_t local_block_beg = locI_extents_[location_id_];
      const uint64_t local_block_end = locI_extents_[location_id_ + 1];
      if (not (gid >= local_block_beg and gid <  local_block_end))
        throw std::logic_error(fname + ": Error mapping index in "
                                       "takers_query_map.");

      takers_local_ids_.push_back(gid - static_cast<int64_t>(local_block_beg));
    }//for gid
  }

  //======================================== Build sendcounts and senddispls
  //Since we now know what the structure will
  //of sending the data we can compute m_sendcounts
  //and m_senddispls, which will be used in a call
  //to MPI_Alltoallv
  sendcounts_.assign(process_count_, 0);
  senddispls_.assign(process_count_, 0);

  uint64_t total_sendcounts = 0;
  for (const auto& pid_vec_pair : taker_query_map)
  {
    const int    pid        = pid_vec_pair.first;
    const size_t num_counts = pid_vec_pair.second.size();
    sendcounts_[pid] = static_cast<int>(num_counts);
    senddispls_[pid] = static_cast<int>(total_sendcounts);
    total_sendcounts += num_counts;
  }
}


//###################################################################
int chi_math::VectorGhostCommunicator::FindOwnerPID(uint64_t global_id) const
{
  if (global_id >= globl_size_)
    throw std::logic_error(
      "chi_math::VectorWithGhosts:FindOwnerPID global_id >= m_globl_size");
  for (int locI=0; locI<static_cast<int>(process_count_); ++locI)
    if (global_id >= locI_extents_[locI] and
        global_id < locI_extents_[locI + 1])
      return locI;
  return int(-1);
}

//###################################################################
void chi_math::VectorGhostCommunicator::
  CommunicateGhostEntries(std::vector<double> &local_vector) const
{
  if (local_vector.size() != (local_size_ + ghost_indices_.size()))
    throw std::logic_error(
      "chi_math::VectorWithGhosts::CommunicateGhostEntries: Vector size "
      "mismatch.");

  //======================================== Build data to be sent
  const size_t amount_to_send = takers_local_ids_.size();
  std::vector<double> send_data;
  send_data.reserve(amount_to_send);
  for (const int64_t local_id : takers_local_ids_)
    send_data.push_back(local_vector[local_id]);

  //======================================== Allocate receive buffer
  const size_t amount_to_recv = ghost_indices_.size();
  std::vector<double> recv_data(amount_to_recv, 0.0);

  //======================================== Communicate
  MPI_Alltoallv(send_data.data(),
                sendcounts_.data(),
                senddispls_.data(),
                MPI_DOUBLE,
                recv_data.data(),
                recvcounts_.data(),
                recvdispls_.data(),
                MPI_DOUBLE,
                comm_);

  //======================================== Populate local vector with ghost
  //                                         data
  for (size_t k=0; k<amount_to_recv; ++k)
  {
    const int64_t gid = ghost_indices_[k];
    const size_t  rid = ghost_ids_to_recv_map_.at(gid);
    local_vector[local_size_ + k] = recv_data[rid];
  }
}


//#########################################################
/**Returns a vector with enough space to hold both `local_size` entries,
 * as specified in the VectorGhostCommunicator-constructor, and entries
 * associated with the ghost-ids. The entries default to zero.*/
std::vector<double> chi_math::VectorGhostCommunicator::
  MakeGhostedVector() const
{
  std::vector<double> vec(local_size_ + ghost_indices_.size(), 0.0);
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
  if (unghosted_vector.size() < local_size_)
    throw std::logic_error(
      "chi_math::VectorGhostCommunicator::MakeGhostedVector: "
      "unghosted_vector.size() < m_local_size");

  std::vector<double> vec(local_size_ + ghost_indices_.size(), 0.0);
  for (size_t i=0; i < local_size_; ++i)
    vec[i] = unghosted_vector[i];
  return vec;
}