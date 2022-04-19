#include "lbs_point_source.h"

#include "chi_log.h"

const chi_mesh::Vector3& lbs::PointSource::Location() const
{
  return m_location;
}

const std::vector<double>& lbs::PointSource::Strength() const
{
  return m_strength;
}

void lbs::PointSource::
  SetOwningCellData(uint64_t owning_cell_local_id,
                    const std::vector<double>& shape_values,
                    const std::vector<double>& node_weights)
{
  m_owning_cell_local_id = owning_cell_local_id;
  m_owning_cell_shape_values = shape_values;
  m_owning_cell_node_weights = node_weights;
  m_owning_cell_set = true;
}

const std::vector<double>& lbs::PointSource::ShapeValues() const
{
  if (not m_owning_cell_set)
    throw std::logic_error("lbs::PointSource::ShapeValues failed because the "
                           "owning cell weights have not been initialized.");

  return m_owning_cell_shape_values;
}

const std::vector<double>& lbs::PointSource::NodeWeights() const
{
  if (not m_owning_cell_set)
    throw std::logic_error("lbs::PointSource::NodeWeights failed because the "
                           "owning cell weights have not been initialized.");

  return m_owning_cell_node_weights;
}


uint64_t lbs::PointSource::OwningCellLocalID() const
{
  if (not m_owning_cell_set)
    throw std::logic_error("lbs::PointSource::OwningCellLocalID failed because the "
                           "owning cell-id has not been initialized.");
  return m_owning_cell_local_id;
}

bool lbs::PointSource::LocallyOwned() const
{
  return m_owning_cell_set;
}

