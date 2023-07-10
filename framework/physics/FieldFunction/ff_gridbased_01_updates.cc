#include "fieldfunction_gridbased.h"

//###################################################################
/**Updates the field data with a STL vector.*/
void chi_physics::FieldFunctionGridBased::
  UpdateFieldVector(const std::vector<double> &field_vector)
{
  if (field_vector.size() < field_vector_.size())
    throw std::logic_error("chi_physics::FieldFunction::UpdateFieldVector: "
                           "Attempted update with a vector of insufficient size.");

  field_vector_ = field_vector;
}

//###################################################################
/**Updates the field data with a PETSc vector.*/
void chi_physics::FieldFunctionGridBased::
  UpdateFieldVector(const Vec& field_vector)
{
  PetscInt n;
  VecGetLocalSize(field_vector, &n);

  if (n < field_vector_.size())
    throw std::logic_error("chi_physics::FieldFunction::UpdateFieldVector: "
                           "Attempted update with a vector of insufficient size.");

  const double* x;
  VecGetArrayRead(field_vector, &x);
  for (size_t i=0; i<n; ++i)
    field_vector_[i] = x[i];
  VecRestoreArrayRead(field_vector, &x);
}