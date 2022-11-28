#include "fieldfunction2.h"

void chi_physics::FieldFunction2::
  UpdateFieldVector(const std::vector<double> &field_vector)
{
  if (field_vector.size() < m_field_vector.size())
    throw std::logic_error("chi_physics::FieldFunction2::UpdateFieldVector: "
                           "Attempted update with vector of different size.");

  m_field_vector = field_vector;
}

const std::vector<double>& chi_physics::FieldFunction2::FieldVector() const
{
  return m_field_vector;
}