#ifndef PARALLEL_VECTOR_H
#define PARALLEL_VECTOR_H

#include "chi_runtime.h"

#include <vector>
#include <cstdint>
#include <map>

#include <mpi.h>


namespace chi_math
{

/**
 * An implementation of a parallel vector.
 */
class ParallelVector
{
public:
  ParallelVector() = delete;

  /**
   * Initialize a parallel vector with the given local and global sizes with
   * the given communicator whose elements are set to zero.
   */
  ParallelVector(const uint64_t local_size,
                 const uint64_t global_size,
                 const MPI_Comm communicator = MPI_COMM_WORLD);

  ParallelVector(const ParallelVector& other)
    : local_size_(other.local_size_),
      global_size_(other.global_size_),
      location_id_(other.location_id_),
      process_count_(other.process_count_),
      values_(other.values_),
      comm_(other.comm_)
  {}

  ParallelVector(ParallelVector&& other)
    : local_size_(other.local_size_),
      global_size_(other.global_size_),
      location_id_(other.location_id_),
      process_count_(other.process_count_),
      values_(other.values_),
      comm_(other.comm_)
  {}

  uint64_t LocalSize() const { return local_size_; }
  virtual uint64_t TotalLocalSize() const { return local_size_; }
  uint64_t GlobalSize() const { return global_size_;}

  /// Return the value of a local entry
  double operator[](const int64_t local_id) const;

  /// Return a reference to a local entry.
  double& operator[](const int64_t local_id);

  /// Set all elements in the parallel vector to the given value.
  virtual void
  Set(const double value) { values_.assign(values_.size(), value); }

  /// Set the parallel vector with the given local vectors.
  virtual void Set(const std::vector<double>& local_vector);

  /**
   * Add the given value to the given global index of the parallel vector.
   *
   * This routine stores the global index-value pair until the Assemble
   * routine is called, when all operations are communicated to the processes
   * that own the corresponding global index.
   */
  void AddValue(const int64_t index, const double value);

  /**
   * Add multiple values to multiple global indices of the parallel vector.
   *
   * This routine simply goes through each global index-value pair and calls
   * AddValue.
   */
  void AddValues(const std::vector<int64_t>& indices,
                 const std::vector<double>& values);

  /**
   * Communicate all operations stored within the operation cache to the
   * corresponding processes that own the respective global indices.
   *
   * This routine clears the operation cache once completed.
   */
  virtual void Assemble();

  /// Print the local vectors to stings.
  std::string PrintStr() const;

protected:
  void DefineParallelStructure();
  int FindOwnerPID(const uint64_t global_id) const;

protected:
  const uint64_t local_size_;
  const uint64_t global_size_;

  int location_id_;
  int process_count_;
  std::vector<uint64_t> extents_;

  std::vector<double> values_;

  using VecAddOp = std::pair<int64_t, double>;
  std::vector<VecAddOp> op_cache_;

  const MPI_Comm comm_;
};

}


#endif //PARALLEL_VECTOR_H
