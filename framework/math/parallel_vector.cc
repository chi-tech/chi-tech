#include "parallel_vector.h"

#include "chi_mpi_utils_map_all2all.h"
#include "data_types/byte_array.h"

#include "chi_log_exceptions.h"

#include <sstream>
#include <stdexcept>


namespace chi_math
{

ParallelVector::ParallelVector(const MPI_Comm communicator)
    : ParallelVector(0, 0, communicator)
{}


ParallelVector::ParallelVector(const uint64_t local_size,
                               const uint64_t global_size,
                               const MPI_Comm communicator)
    : communicator_(communicator)
{
  // Get the processor ID and the number of processors
  MPI_Comm_rank(communicator_, &location_id_);
  MPI_Comm_size(communicator_, &process_count_);

  // Reinitialize the vector
  Reinit(local_size, global_size);
}


void ParallelVector::Clear()
{
  local_size_ = 0;
  global_size_ = 0;
  extents_.clear();
  values_.clear();
  op_cache_.clear();
}


void ParallelVector::Reinit(const uint64_t local_size,
                            const uint64_t global_size)
{
  Clear();

  local_size_ = local_size;
  global_size_ = global_size;

  // Set the local vector to zero
  values_.assign(local_size_, 0.0);

  // Get the local vector sizes per processor
  std::vector<uint64_t> locJ_local_size(process_count_, 0);
  MPI_Allgather(&local_size,            //sendbuf
                1, MPI_UINT64_T,        //sendcount + sendtype
                locJ_local_size.data(), //recvbuf
                1, MPI_UINT64_T,        //recvcount + recvtype
                communicator_);         //communicator

  // With the vector sizes per processor, now the offsets for each
  // processor can be defined using a cumulative sum per processor.
  // This allows for the determination of whether a global index is
  // locally owned or not.
  extents_.assign(process_count_ + 1, 0);
  for (size_t locJ = 1; locJ < process_count_; ++locJ)
    extents_[locJ] = extents_[locJ - 1] + locJ_local_size[locJ - 1];
  extents_[process_count_] =
      extents_[process_count_ - 1] + locJ_local_size.back();

}


double ParallelVector::operator[](const int64_t i) const
{
  ChiInvalidArgumentIf(
      i < 0 or i >= local_size_,
      std::string(__FUNCTION__) + ": " +
      "Invalid local ID encountered on process " +
      std::to_string(location_id_) + ". Local IDs must be in the range " +
      "[0, local_size_ - 1].");

  return values_[i];
}


double& ParallelVector::operator[](const int64_t i)
{
  ChiInvalidArgumentIf(
      i < 0 or i >= local_size_,
      std::string(__FUNCTION__) + ": " +
      "Invalid local ID encountered on process " +
      std::to_string(location_id_) + ". Local IDs must be in the range " +
      "[0, local_size_ - 1].");

  return values_[i];
}


void ParallelVector::Set(const double value)
{
  values_.assign(local_size_, value);
}


void ParallelVector::AddValue(const int64_t index, const double value)
{
  ChiInvalidArgumentIf(
      index < 0 or index >= global_size_,
      std::string(__FUNCTION__) + ": " +
      "Invalid global index encountered. Global indices cannot be "
      "negative or exceed the global vector size. Provided index was " +
      std::to_string(index) + " and global size is " +
      std::to_string(global_size_) + ".");

  op_cache_.push_back({index, value});
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
    const int pid = FindProcessID(global_id);
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


int ParallelVector::FindProcessID(const uint64_t global_id) const
{
  ChiInvalidArgumentIf(
      global_id >= global_size_,
      std::string(__FUNCTION__) + ": " +
      "The specified global ID is greater than the global vector size.");

  for (int locJ = 0; locJ < process_count_; ++locJ)
    if (global_id >= extents_[locJ] and
        global_id < extents_[locJ + 1])
      return locJ;
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
