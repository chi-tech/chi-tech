#ifndef PARALLEL_VECTOR_H
#define PARALLEL_VECTOR_H

#include "ParallelVector.h"
#include "math/chi_math.h"

#include <vector>
#include <cstdint>
#include <map>

#include <mpi.h>

namespace chi_math
{

/**
 * An implementation of a parallel vector.
 */
class ParallelSTLVector : public ParallelVector
{
public:
  /**
   * Initialize a parallel vector with the given local and global sizes with
   * the given communicator whose entries are set to zero.
   */
  ParallelSTLVector(uint64_t local_size,
                 uint64_t global_size,
                 MPI_Comm communicator);

  /**Copy constructor.*/
  ParallelSTLVector(const ParallelSTLVector& other);

  /**Move constructor.*/
  ParallelSTLVector(ParallelSTLVector&& other) noexcept;

  std::unique_ptr<ParallelVector> MakeCopy() const override;
  std::unique_ptr<ParallelVector> MakeNewVector() const override;

  /** Returns the raw stl-vector associated with the local data of this
   * vector.*/
  std::vector<double>& RawValues() override { return values_; }

  /**
   * Read only accessor to the entry at the given local index of
   * the local vector.
   *
   * \note This accessor allows access to all locally stored elements,
   *       including any data beyond local_size_ that may exist in derived
   *       classes.
   */
  double operator[](int64_t local_id) const override;

  /**
   * Read/write accessor to the entry at the given local index of
   * the local vector.
   *
   * \note This accessor allows access to all locally stored elements,
   *       including any data beyond local_size_ that may exist in derived
   *       classes.
   */
  double& operator[](int64_t local_id) override;

  /// Return a vector containing the locally owned data.
  std::vector<double> MakeLocalVector() override;

  /**
   * Set the entries of the locally owned portion of the parallel vector
   * to the given value.
   */
  void Set(const double value) override
  {
    values_.assign(values_.size(), value);
  }

  /**
   * Set the entries of the locally owned portion of the parallel vector
   * to the given STL vector.
   */
  void Set(const std::vector<double>& local_vector) override;

  /**
   * Define a set or add operation for the given global id-value pair
   *
   * This routine adds the global id-value pair to the set operation cache,
   * which upon execution of Assemble, communicates the operations to the
   * appropriate process.
   */
  void SetValue(int64_t global_id, double value, VecOpType op_type) override;

  /**
   * Group multiple operations into a single call.
   *
   * This routine goes through the given global id-value pairs and calls
   * SetValue for each.
   */
  void SetValues(const std::vector<int64_t>& global_ids,
                 const std::vector<double>& values,
                 VecOpType op_type) override;

  /**In place adding of vectors. The sizes must be compatible.*/
  void operator+=(const ParallelVector& y) override;

  /**In place adding of vectors. The sizes must be compatible.*/
  void operator+=(const ParallelSTLVector& y);

  /**Adds a vector multiplied by scalar a. Optimized for a=1.0 and -1.0*/
  void PlusAY(double a, const ParallelVector& y) override;

  /**Sets the local values of one vector equal to another. The sizes must be
   * compatible.*/
  void CopyValues(const ParallelVector& y) override;

  /**Sets the local values of one vector equal to another. The sizes must be
   * compatible.*/
  void CopyValues(const ParallelSTLVector& y);

  void BlockCopyLocalValues(const ParallelVector& y,
                            int64_t y_offset,
                            int64_t local_offset,
                            int64_t num_values) override;

  /**Returns the specified norm of the vector.*/
  double ComputeNorm(chi_math::NormType norm_type) const override;

  /**
   * Communicate all operations stored within the operation cache to the
   * corresponding processes that own the respective global indices, and
   * apply the operations.
   *
   * This routine clears the respective operation cache once completed.
   */
  void Assemble() override;

  /// Print the local vectors to stings.
  std::string PrintStr() const override;

protected:
  const std::vector<uint64_t> extents_;

  std::vector<double> values_;

  using Operation = std::pair<int64_t, double>;
  std::vector<Operation> set_cache_;
  std::vector<Operation> add_cache_;

private:
  static std::vector<uint64_t>
  DefineExtents(uint64_t local_size, int comm_size, MPI_Comm communicator);
  int FindOwnerPID(uint64_t global_id) const;
};

} // namespace chi_math

#endif // PARALLEL_VECTOR_H
