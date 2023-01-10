#ifndef CHITECH_LBS_POINT_SOURCE_H
#define CHITECH_LBS_POINT_SOURCE_H

#include <utility>

#include "ChiMesh/chi_mesh.h"

namespace lbs
{

class PointSource
{
private:
  const chi_mesh::Vector3 m_location;
  const std::vector<double> m_groupwise_strength;

public:
  struct ContainingCellInfo
  {
    double volume_weight;
    uint64_t cell_local_id;
    std::vector<double> shape_values;
    std::vector<double> node_weights;
  };
private:
  std::vector<ContainingCellInfo> m_containing_cells;

public:
  PointSource(const chi_mesh::Vector3& location,
              std::vector<double>  strength) :
              m_location(location),
              m_groupwise_strength(std::move(strength))
  {}

  const chi_mesh::Vector3& Location() const;

  const std::vector<double>& Strength() const;

  void AddContainingCellInfo(double volume_weight,
                             uint64_t cell_local_id,
                             std::vector<double> shape_values,
                             std::vector<double> node_weights);

  const std::vector<ContainingCellInfo>& ContainingCellsInfo() const;

  void ClearInitializedInfo();
};

}//namespace lbs

#endif //CHITECH_LBS_POINT_SOURCE_H
