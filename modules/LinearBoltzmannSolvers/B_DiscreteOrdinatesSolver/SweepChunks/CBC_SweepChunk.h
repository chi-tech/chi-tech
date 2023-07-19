#ifndef CHITECH_CBC_SWEEPCHUNK_H
#define CHITECH_CBC_SWEEPCHUNK_H

#include "SweepChunk.h"

namespace lbs
{

class CBC_FLUDS;
class CBC_ASynchronousCommunicator;

struct CBC_SweepDependencyInterface : public SweepDependencyInterface
{
  CBC_FLUDS* fluds_ = nullptr;
  CBC_ASynchronousCommunicator* communicator_;

  // Set using SetupIncomingFace
  const chi_mesh::sweep_management::FaceNodalMapping* face_nodal_mapping_ =
    nullptr;

  const double* GetUpwindPsi(int face_node_local_idx) const override;
  double* GetDownwindPsi(int face_node_local_idx) const override;
  void SetupIncomingFace(int face_id,
                         size_t num_face_nodes,
                         uint64_t neighbor_id,
                         bool on_local_face,
                         bool on_boundary) override;
  void SetupOutgoingFace(int face_id,
                         size_t num_face_nodes,
                         uint64_t neighbor_id,
                         bool on_local_face,
                         bool on_boundary,
                         int locality) override;
};

class CBC_SweepChunk : public SweepChunk
{
public:
  CBC_SweepChunk(std::vector<double>& destination_phi,
                 std::vector<double>& destination_psi,
                 const chi_mesh::MeshContinuum& grid,
                 const chi_math::SpatialDiscretization& discretization,
                 const std::vector<UnitCellMatrices>& unit_cell_matrices,
                 std::vector<lbs::CellLBSView>& cell_transport_views,
                 const std::vector<double>& source_moments,
                 const LBSGroupset& groupset,
                 const std::map<int, XSPtr>& xs,
                 int num_moments,
                 int max_num_cell_dofs);

  void SetCell(chi_mesh::Cell const* cell_ptr) override;

  void Sweep(chi_mesh::sweep_management::AngleSet& angle_set) override;

protected:
  chi_mesh::Cell const* cell_ptr_ = nullptr;
  uint64_t cell_local_id_ = 0;
};

} // namespace lbs

#endif // CHITECH_CBC_SWEEPCHUNK_H
