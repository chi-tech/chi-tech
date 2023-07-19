#include "CBC_FLUDS.h"

#include "math/SpatialDiscretization/spatial_discretization.h"
#include "mesh/SweepUtilities/SPDS/SPDS.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

namespace lbs
{

CBC_FLUDS::CBC_FLUDS(size_t num_groups,
                     size_t num_angles,
                     const CBC_FLUDSCommonData& common_data,
                     std::vector<double>& local_psi_data,
                     const chi_math::UnknownManager& psi_uk_man,
                     const chi_math::SpatialDiscretization& sdm)
  : chi_mesh::sweep_management::FLUDS(
      num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    local_psi_data_(local_psi_data),
    psi_uk_man_(psi_uk_man),
    sdm_(sdm)
{
}

const chi_mesh::sweep_management::FLUDSCommonData& CBC_FLUDS::CommonData() const
{
  return common_data_;
}

const double* CBC_FLUDS::GetUpwindPsi(uint64_t cell_global_id,
                                      unsigned int cell_node,
                                      unsigned int angle_id,
                                      size_t g)
{
  const auto& cell = spds_.Grid().cells[cell_global_id];
  const auto dof_map =
    sdm_.MapDOFLocal(cell, cell_node, psi_uk_man_, angle_id, g);

  return &local_psi_data_.get()[dof_map];
}

const double* CBC_FLUDS::GetNLUpwindPsi(uint64_t cell_global_id,
                                        unsigned int face_id,
                                        unsigned int face_node_mapped,
                                        unsigned int angle_set_index,
                                        size_t group_angle_stride,
                                        size_t group_stride)
{
  const size_t dof_map =
    face_node_mapped * group_angle_stride + angle_set_index * group_stride;

  return &deplocs_outgoing_messages_.at({cell_global_id, face_id}).at(dof_map);
}

} // namespace lbs