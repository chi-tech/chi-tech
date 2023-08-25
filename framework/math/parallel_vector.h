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
  /// Initialize an empty parallel vector.
  ParallelVector(const MPI_Comm communicator = MPI_COMM_WORLD);

  /// Initialize a zero parallel vector with the given local and global sizes.
  ParallelVector(const uint64_t local_size,
                 const uint64_t global_size,
                 const MPI_Comm communicator = MPI_COMM_WORLD);

  /// Return the local size of the vector.
  uint64_t LocalSize() const { return local_size_; }

  /// Return the global size of the vector.
  uint64_t GlobalSize() const { return global_size_;}

  /// Clear the parallel vector, returning it to an uninitialized state.
  virtual void Clear();

  /**
   * Reinitialize the parallel vector to the given local and global sizes.
   *
   * This routine clears all data in the current parallel vectors and
   * results in a zero vector.
   */
  void Reinit(const uint64_t local_size,
              const uint64_t global_size);

  /**
   * Set all elements in the parallel vector to the given value.
   */
  void Set(const double value);

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

  /// Return the <tt>i</tt>'th element of the local vector.
  double operator[](const int64_t i) const;

  /// Return a reference to the <tt>i</tt>'th element of the local vector.
  double& operator[](const int64_t i);

  std::string PrintStr() const;

protected:
  int FindProcessID(const uint64_t global_id) const;

protected:
  uint64_t local_size_;
  uint64_t global_size_;
  const MPI_Comm communicator_;

  int location_id_ = 0;
  int process_count_ = 0;
  std::vector<uint64_t> extents_;

  std::vector<double> values_;

  using VecAddOp = std::pair<int64_t, double>;
  std::vector<VecAddOp> op_cache_;

};

}


#endif //PARALLEL_VECTOR_H
