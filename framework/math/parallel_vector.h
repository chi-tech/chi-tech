#ifndef PARALLEL_VECTOR_H
#define PARALLEL_VECTOR_H

#include "chi_runtime.h"

#include <vector>
#include <cstdint>
#include <map>

#include <mpi.h>


namespace chi_math
{

enum class OperationType : short
{
  SET_VALUE = 1,
  ADD_VALUE = 2
};


/**
 * An implementation of a parallel vector.
 */
class ParallelVector
{
public:
  using iterator = std::vector<double>::iterator;
  using const_iterator = std::vector<double>::const_iterator;

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

  double operator[](const int64_t local_id) const;
  double& operator[](const int64_t local_id);

  iterator begin() { return values_.begin(); }
  iterator end() { return values_.begin() + local_size_; }

  const_iterator begin() const { return values_.begin(); }
  const_iterator end() const { return values_.end(); }

  std::vector<double> MakeLocalVector();

  /// Set all elements of the local parallel vector to the given value.
  virtual void
  Set(const double value) { values_.assign(values_.size(), value); }

  /// Set all elements of the local parallel vector with an STL vector.
  virtual void Set(const std::vector<double>& local_vector);

  /**
   * Define a set or add operation for the given global id-value pair
   *
   * This routine adds the global id-value pair to the set operation cache,
   * which upon execution of Assemble, communicates the operations to the
   * appropriate process.
   */
  void SetValue(const int64_t global_id,
                const double value,
                const OperationType op_type);

  /**
   * Group multiple operations into a single call.
   *
   * This routine goes through the given global id-value pairs and calls
   * SetValue for each.
   */
  void SetValues(const std::vector<int64_t>& global_ids,
                 const std::vector<double>& values,
                 const OperationType op_type);

  /**
   * Communicate all operations stored within the operation cache to the
   * corresponding processes that own the respective global indices, and
   * apply the operations.
   *
   * This routine clears the respective operation cache once completed.
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

  using Operation = std::pair<int64_t, double>;
  std::vector<Operation> set_cache_;
  std::vector<Operation> add_cache_;

  const MPI_Comm comm_;
};

}


#endif //PARALLEL_VECTOR_H
