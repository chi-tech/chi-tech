#include "math/SpatialDiscretization/CellMappings/FiniteVolumeMapping.h"

#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

namespace chi_math::cell_mapping
{

finite_element::InternalQuadraturePointData
FiniteVolumeMapping::MakeInternalQuadraturePointData() const
{
  return finite_element::InternalQuadraturePointData(
    {0},
    {{cell_.centroid_}},
    {{1.0}},
    {{chi_mesh::Vector3(0, 0, 0)}},
    {volume_},
    face_node_mappings_,
    num_nodes_);
}

finite_element::FaceQuadraturePointData
FiniteVolumeMapping::MakeFaceQuadraturePointData(size_t face_index) const
{
  return finite_element::FaceQuadraturePointData({0},
                                                 {{chi_mesh::Vector3(0, 0, 0)}},
                                                 {{1.0}},
                                                 {{chi_mesh::Vector3(0, 0, 0)}},
                                                 {areas_[face_index]},
                                                 {{chi_mesh::Vector3(0, 0, 0)}},
                                                 face_node_mappings_,
                                                 1);
}

} // namespace chi_math
