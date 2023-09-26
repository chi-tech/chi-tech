#ifndef CHITECH_PARALLELVECTOR_H
#define CHITECH_PARALLELVECTOR_H

#include "math/chi_math.h"

#include <cstdint>

namespace chi_math
{

enum class VecOpType : short
{
  SET_VALUE = 1,
  ADD_VALUE = 2
};

/**
 * An abstract implementation of a parallel vector.
 */
class ParallelVector
{
public:
  /**
   * Initialize a parallel vector with the given local and global sizes with
   * the given communicator whose entries are set to zero.
   */
  ParallelVector(uint64_t local_size,
                 uint64_t global_size,
                 MPI_Comm communicator);

  /**Copy constructor.*/
  ParallelVector(const ParallelVector& other);

  /**Move constructor.*/
  ParallelVector(ParallelVector&& other) noexcept;

  /**Creates a copy of the vector datastructures AND values.
  * This routine requires no communication.*/
  virtual std::unique_ptr<ParallelVector> MakeCopy() const = 0;

  /**Creates a copy of the vector datastructures but NOT the values. The
  * values are defaulted to zero. This routine requires no communication.*/
  virtual std::unique_ptr<ParallelVector> MakeNewVector() const = 0;

  /** Returns the raw stl-vector associated with the local data of this
   * vector.*/
  virtual std::vector<double>& RawValues() = 0;

  /// Return the size of the locally owned portion of the parallel vector.
  uint64_t LocalSize() const { return local_size_; }

  /// Return the global size of the parallel vector.
  uint64_t GlobalSize() const { return global_size_; }

  /**
   * Read only accessor to the entry at the given local index of
   * the local vector.
   *
   * \note This accessor allows access to all locally stored elements,
   *       including any data beyond local_size_ that may exist in derived
   *       classes.
   */
  virtual double operator[](int64_t local_id) const = 0;

  /**
   * Read/write accessor to the entry at the given local index of
   * the local vector.
   *
   * \note This accessor allows access to all locally stored elements,
   *       including any data beyond local_size_ that may exist in derived
   *       classes.
   */
  virtual double& operator[](int64_t local_id) = 0;

  /// Return a vector containing the locally owned data.
  virtual std::vector<double> MakeLocalVector() = 0;

  /**
   * Set the entries of the locally owned portion of the parallel vector
   * to the given value.
   */
  virtual void Set(double value) = 0;

  /**
   * Set the entries of the locally owned portion of the parallel vector
   * to the given STL vector.
   */
  virtual void Set(const std::vector<double>& local_vector) = 0;

  /**
   * Define a set or add operation for the given global id-value pair
   *
   * This routine adds the global id-value pair to the set operation cache,
   * which upon execution of Assemble, communicates the operations to the
   * appropriate process.
   */
  virtual void SetValue(int64_t global_id, double value, VecOpType op_type) = 0;

  /**
   * Group multiple operations into a single call.
   *
   * This routine goes through the given global id-value pairs and calls
   * SetValue for each.
   */
  virtual void SetValues(const std::vector<int64_t>& global_ids,
                         const std::vector<double>& values,
                         VecOpType op_type) = 0;

  /**In place adding of vectors. The sizes must be compatible.*/
  virtual void operator+=(const ParallelVector& y) = 0;

  /**Adds a vector multiplied by scalar a. Optimized for a=1.0 and -1.0*/
  virtual void PlusAY(double a, const ParallelVector& y) = 0;

  /**Sets the local values of one vector equal to another. The sizes must be
   * compatible.*/
  virtual void CopyValues(const ParallelVector& y) = 0;

  /**Copies a contiguous block of local data (num_values entries) from the
  * source vector (starting at y_offset) to the
   * current vector starting at local_offset. */
  virtual void BlockCopyLocalValues(const ParallelVector& y,
                                    int64_t y_offset,
                                    int64_t local_offset,
                                    int64_t num_values) = 0;

  /**Copies a contiguous block of data from the source STL vector to the
   * current vector starting at local_offset. The input STL vector must have
   * exactly num_values entries.
   * */
  virtual void BlockCopyLocalValues(const std::vector<double>& y,
                                    int64_t local_offset,
                                    int64_t num_values) = 0;

  /**Returns the specified norm of the vector.*/
  virtual double ComputeNorm(chi_math::NormType norm_type) const = 0;

  /**
   * Communicate all operations stored within the operation cache to the
   * corresponding processes that own the respective global indices, and
   * apply the operations.
   *
   * This routine clears the respective operation cache once completed.
   */
  virtual void Assemble() = 0;

  /// Print the local vectors to stings.
  virtual std::string PrintStr() const = 0;

  /**
   * Communicate the current ghost entries, if applicable, to all other
   * processes to update the locally stored ghost data.
   */
  virtual void CommunicateGhostEntries(){};

  virtual ~ParallelVector() = default;

protected:
  const uint64_t local_size_;
  const uint64_t global_size_;

  const int location_id_;
  const int process_count_;
  const MPI_Comm comm_;
};

} // namespace chi_math

#endif // CHITECH_PARALLELVECTOR_H
