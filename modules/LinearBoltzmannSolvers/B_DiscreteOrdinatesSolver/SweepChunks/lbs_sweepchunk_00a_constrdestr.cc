#include "lbs_sweepchunk.h"

namespace lbs
{

LBSSweepChunk::LBSSweepChunk(
  const chi_mesh::MeshContinuum& grid,
  const chi_math::SpatialDiscretization& discretization,
  const std::vector<UnitCellMatrices>& unit_cell_matrices,
  std::vector<lbs::CellLBSView>& cell_transport_views,
  std::vector<double>& destination_phi,
  std::vector<double>& destination_psi,
  const std::vector<double>& source_moments,
  LBSGroupset& groupset,
  const std::map<int, XSPtr>& xs,
  int num_moments,
  int max_num_cell_dofs)
  : SweepChunk(destination_phi, destination_psi),
    grid_(grid),
    grid_fe_view_(discretization),
    unit_cell_matrices_(unit_cell_matrices),
    grid_transport_view_(cell_transport_views),
    q_moments_(source_moments),
    groupset_(groupset),
    xs_(xs),
    num_moments_(num_moments),
    num_groups_(groupset.groups_.size()),
    save_angular_flux_(!destination_psi.empty())
{
  Amat_.resize(max_num_cell_dofs, std::vector<double>(max_num_cell_dofs));
  Atemp_.resize(max_num_cell_dofs, std::vector<double>(max_num_cell_dofs));
  b_.resize(num_groups_, std::vector<double>(max_num_cell_dofs, 0.0));
  source_.resize(max_num_cell_dofs, 0.0);
}

} // namespace lbs