#include "lbs_point_source.h"

#include "chi_log.h"

/**Returns the location of the source.*/
const chi_mesh::Vector3& lbs::PointSource::Location() const
{
  return m_location;
}

/**Returns the groupwise strength of the source.*/
const std::vector<double>& lbs::PointSource::Strength() const
{
  return m_groupwise_strength;
}

/**Adds the given information as `ContainingCellInfo` structure.
 * This info is necessary to properly normalize a source.*/
void lbs::PointSource::
  AddContainingCellInfo(double volume_weight,
                        uint64_t cell_local_id,
                        std::vector<double> shape_values,
                        std::vector<double> node_weights)
{
  m_containing_cells.push_back(ContainingCellInfo{volume_weight,
                                                  cell_local_id,
                                                  std::move(shape_values),
                                                  std::move(node_weights)});
}

/**Get a list of cell-information for cells containing/sharing this source.*/
const std::vector<lbs::PointSource::ContainingCellInfo>&
  lbs::PointSource::ContainingCellsInfo() const
{
  return m_containing_cells;
}

/**Clears all initialization info so that the source can be reinitialized.*/
void lbs::PointSource::ClearInitializedInfo()
{
  m_containing_cells.clear();
}