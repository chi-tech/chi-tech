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
  /**
   * Initialize a parallel vector with the given local and global sizes with
   * the given communicator whose entries are set to zero.
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

  /// Return the size of the locally owned portion of the parallel vector.
  uint64_t LocalSize() const { return local_size_; }

  /// Return the global size of the parallel vector.
  uint64_t GlobalSize() const { return global_size_;}

  /**
   * Read only accessor to the entry at the given local index of
   * the local vector.
   *
   * \note This accessor allows access to all locally stored elements,
   *       including any data beyond local_size_ that may exist in derived
   *       classes.
   */
  double operator[](const int64_t local_id) const;

  /**
   * Read/write accessor to the entry at the given local index of
   * the local vector.
   *
   * \note This accessor only allows access to the locally owned entries
   *       of the parallel vector, and does not allow access for any data
   *       beyond local_size_ that may exist in derived classes.
   */
  double& operator[](const int64_t local_id);

  /// Return a vector containing the locally owned data.
  std::vector<double> MakeLocalVector();

  /**
   * Set the entries of the locally owned portion of the parallel vector
   * to the given value.
   */
  void Set(const double value) { values_.assign(values_.size(), value); }

  /**
   * Set the entries of the locally owned portion of the parallel vector
   * to the given STL vector.
   */
  void Set(const std::vector<double>& local_vector);

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
