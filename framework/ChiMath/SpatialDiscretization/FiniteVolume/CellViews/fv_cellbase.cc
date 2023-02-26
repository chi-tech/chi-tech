#include "fv_cellbase.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

void chi_math::CellFVValues::InitializeVolumeQuadraturePointData(
  finite_element::InternalQuadraturePointData &internal_data) const
{
  internal_data.InitializeData({0},
                               {{cell_.centroid}},
                               {{1.0}},
                               {{chi_mesh::Vector3(0, 0, 0)}},
                               {volume_},
                               face_node_mappings_,
                               num_nodes_);
}

void chi_math::CellFVValues::
  InitializeFaceQuadraturePointData(
    unsigned int face,
    finite_element::FaceQuadraturePointData &faces_qp_data) const
{
  faces_qp_data.InitializeData({0},
                               {{chi_mesh::Vector3(0, 0, 0)}},
                               {{1.0}},
                               {{chi_mesh::Vector3(0, 0, 0)}},
                               {areas_[face]},
                               {{chi_mesh::Vector3(0, 0, 0)}},
                               face_node_mappings_,
                               1);
}