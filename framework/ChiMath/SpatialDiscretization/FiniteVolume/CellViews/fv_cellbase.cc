#include "fv_cellbase.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

void chi_math::CellFVValues::InitializeVolumeQuadraturePointData(
  finite_element::InternalQuadraturePointData &internal_data) const
{
  internal_data.InitializeData({0},
                               {{m_cell.centroid_}},
                               {{1.0}},
                               {{chi_mesh::Vector3(0, 0, 0)}},
                               {m_volume},
                               face_node_mappings,
                               m_num_nodes);
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
                               {m_areas[face]},
                               {{chi_mesh::Vector3(0, 0, 0)}},
                               face_node_mappings,
                               1);
}