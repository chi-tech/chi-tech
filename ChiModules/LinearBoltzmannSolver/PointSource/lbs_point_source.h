#ifndef CHITECH_LBS_POINT_SOURCE_H
#define CHITECH_LBS_POINT_SOURCE_H

#include "ChiMesh/chi_mesh.h"

namespace lbs
{

class PointSource
{
private:
  const chi_mesh::Vector3 m_location;
  const std::vector<double> m_strength;

  bool m_owning_cell_set = false;
  uint64_t m_owning_cell_local_id = 0;

  std::vector<double> m_owning_cell_shape_values;
  std::vector<double> m_owning_cell_node_weights;

public:
  PointSource(const chi_mesh::Vector3& location,
              const std::vector<double>& strength) :
              m_location(location),
              m_strength(strength)
  {}

  const chi_mesh::Vector3& Location() const;

  const std::vector<double>& Strength() const;

  void SetOwningCellData(uint64_t owning_cell_local_id,
                         const std::vector<double>& shape_values,
                         const std::vector<double>& node_weights);

  const std::vector<double>& ShapeValues() const;
  const std::vector<double>& NodeWeights() const;

  uint64_t OwningCellLocalID() const;

  bool LocallyOwned() const;
};

}//namespace lbs

#endif //CHITECH_LBS_POINT_SOURCE_H
