#include "parallel_vector.h"

#include "chi_mpi_utils_map_all2all.h"
#include "data_types/byte_array.h"

#include "chi_log_exceptions.h"

#include <sstream>
#include <stdexcept>


namespace chi_math
{

ParallelVector::ParallelVector(const uint64_t local_size,
                               const uint64_t global_size,
                               const MPI_Comm communicator)
  : local_size_(local_size),
    global_size_(global_size),
    comm_(communicator)
{
  // Get the processor ID and the number of processors
  MPI_Comm_rank(comm_, &location_id_);
  MPI_Comm_size(comm_, &process_count_);

  DefineParallelStructure();
  values_.assign(local_size_, 0.0);
}


std::vector<double> ParallelVector::MakeLocalVector()
{
  return std::vector<double>(values_.begin(),
                             values_.begin() + local_size_);
}


double ParallelVector::operator[](const int64_t local_id) const
{
  ChiInvalidArgumentIf(
      local_id < 0 or local_id >= local_size_,
      std::string(__FUNCTION__) + ": Invalid local index provided.");
  return values_[local_id];
}


double& ParallelVector::operator[](const int64_t local_id)
{
  ChiInvalidArgumentIf(
      local_id < 0 or local_id >= local_size_,
      std::string(__FUNCTION__) + ": Invalid local index provided.");
  return values_[local_id];
}


void ParallelVector::Set(const std::vector<double>& local_vector)
{
  ChiInvalidArgumentIf(
      local_vector.size() != local_size_,
      std::string(__FUNCTION__) + ": Incompatible local vector size.");
  values_ = local_vector;
}


void ParallelVector::AddValue(const int64_t index, const double value)
{
  if (index >= extents_[location_id_] and
      index < extents_[location_id_ + 1])
    values_[index - extents_[location_id_]] += value;
  else if (index >= 0 and index < global_size_)
    op_cache_.push_back({index, value});
  else
    ChiInvalidArgument("Invalid global index encountered. Global indices "
                       "must be in the range [0, this->GlobalSize()].");
}


void ParallelVector::AddValues(const std::vector<int64_t>& indices,
                               const std::vector<double>& values)
{
  ChiInvalidArgumentIf(
      indices.size() != values.size(),
      std::string(__FUNCTION__) + ": " +
      "Size mismatch between indices and values.");

  for (size_t i = 0; i < indices.size(); ++i)
    AddValue(indices[i], values[i]);
}


void ParallelVector::Assemble()
{
  // Define a mapping between processes and serialized operations that
  // need to be communicated.
  std::map<int, chi_data_types::ByteArray> pid_send_map;
  for (const auto& [global_id, value] : op_cache_)
  {
    const int pid = FindOwnerPID(global_id);
    auto& byte_array = pid_send_map[pid];
    byte_array.Write(global_id);
    byte_array.Write(value);
  }

  // Convert this mapping to a vector of bytes.
  std::map<int, std::vector<std::byte>> pid_send_map_bytes;
  for (const auto& [pid, byte_array] : pid_send_map)
    pid_send_map_bytes[pid] = byte_array.Data();

  // Next, communicate what is going where to all other processes so that
  // each process knows what it will be receiving.
  std::map<int, std::vector<std::byte>> pid_recv_map_bytes =
      chi_mpi_utils::MapAllToAll(pid_send_map_bytes, MPI_BYTE);

  // Now, deserialize the received data and contribute it to the vector
  for (const auto& [pid, byte_vector] : pid_recv_map_bytes)
  {
    const auto packet_size = sizeof(int64_t) + sizeof(double);

    ChiLogicalErrorIf(
        byte_vector.size() % packet_size != 0,
        std::string(__FUNCTION__) + ": " +
        "Unrecognized received operations. Operations are serialized with "
        "an int64_t and double, but the received packet from process " +
        std::to_string(pid) + " on process " + std::to_string(location_id_) +
        " is not an integer multiple of the size of an int64_t and double.");

    const size_t num_ops = byte_vector.size() / packet_size ;
    chi_data_types::ByteArray byte_array(byte_vector);
    for (size_t k = 0; k < num_ops; ++k)
    {
      const int64_t global_id = byte_array.Read<int64_t>();
      const double value = byte_array.Read<double>();

      // Check that the global ID is in fact valid for this process
      const int64_t local_id = global_id - extents_[location_id_];

      ChiLogicalErrorIf(
          local_id < 0 or local_id >= local_size_,
          std::string(__FUNCTION__) + ": " +
          "A non-local global ID was received by process " +
          std::to_string(location_id_) + " by process " +
          std::to_string(pid) + " during vector assembly.");

      // Contribute to the local vector
      values_[local_id] += value;
    }
  }

  // Finally, clear the operation cache
  op_cache_.clear();
}


void ParallelVector::DefineParallelStructure()
{
  // Get the local vector sizes per processor
  std::vector<uint64_t> local_sizes(process_count_, 0);
  MPI_Allgather(&local_size_,            //sendbuf
                1, MPI_UINT64_T,        //sendcount + sendtype
                local_sizes.data(), //recvbuf
                1, MPI_UINT64_T,        //recvcount + recvtype
                comm_);         //communicator

  // With the vector sizes per processor, now the offsets for each
  // processor can be defined using a cumulative sum per processor.
  // This allows for the determination of whether a global index is
  // locally owned or not.
  extents_.assign(process_count_ + 1, 0);
  for (size_t p = 1; p < process_count_; ++p)
    extents_[p] = extents_[p - 1] + local_sizes[p - 1];
  extents_[process_count_] =
      extents_[process_count_ - 1] + local_sizes.back();
}


int ParallelVector::FindOwnerPID(const uint64_t global_id) const
{
  ChiInvalidArgumentIf(
      global_id >= global_size_,
      std::string(__FUNCTION__) + ": " +
      "The specified global ID is greater than the global vector size.");

  for (int p = 0; p < process_count_; ++p)
    if (global_id >= extents_[p] and
        global_id < extents_[p + 1])
      return p;
  return -1;

}


std::string ParallelVector::PrintStr() const
{
  std::stringstream ss;
  for (size_t i = 0; i < values_.size(); ++i)
    ss << (i == 0 ? "[" : "") << values_[i]
       << (i < values_.size() - 1 ? " " : "]");
  return ss.str();
}

}
