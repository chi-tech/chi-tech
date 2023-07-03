#ifndef CHI_SPDS_H
#define CHI_SPDS_H

#include "mesh/SweepUtilities/SPLS/SPLS.h"

#include <memory>

namespace chi_mesh::sweep_management
{

// ###################################################################
/**Contains multiple levels*/
class SPDS
{
public:
  SPDS(const chi_mesh::Vector3& in_omega,
       const chi_mesh::MeshContinuum& in_grid) : omega_(in_omega), grid_(in_grid)
  {
  }

  const chi_mesh::MeshContinuum& Grid() const { return grid_; }
  const chi_mesh::Vector3& Omega() const { return omega_; }
  const SPLS& GetSPLS() const { return spls_; }
  const std::vector<STDG>& GetGlobalSweepPlanes() const
  {
    return global_sweep_planes_;
  }
  typedef std::vector<int> VecInt;
  const VecInt& GetLocationDependencies() const
  {
    return location_dependencies_;
  }
  const VecInt& GetLocationSuccessors() const { return location_successors_; }
  const VecInt& GetDelayedLocationDependencies() const
  {
    return delayed_location_dependencies_;
  }
  const VecInt& GetDelayedLocationSuccessors() const
  {
    return delayed_location_successors_;
  }
  const std::vector<std::pair<int, int>>& GetLocalCyclicDependencies() const
  {
    return local_cyclic_dependencies_;
  }
  const std::vector<std::vector<FaceOrientation>>& CellFaceOrientations() const
  {
    return cell_face_orientations_;
  }

  int MapLocJToPrelocI(int locJ) const;
  int MapLocJToDeplocI(int locJ) const;

protected:
  chi_mesh::Vector3 omega_;

  const chi_mesh::MeshContinuum& grid_;

  SPLS spls_;
  std::vector<STDG> global_sweep_planes_; ///< Processor sweep planes
  std::vector<int> location_dependencies_;
  std::vector<int> location_successors_;
  std::vector<int> delayed_location_dependencies_;
  std::vector<int> delayed_location_successors_;

  std::vector<std::pair<int, int>> local_cyclic_dependencies_;

  std::vector<std::vector<FaceOrientation>> cell_face_orientations_;

  void BuildTaskDependencyGraph(
    const std::vector<std::vector<int>>& global_dependencies,
    bool cycle_allowance_flag);
};

} // namespace chi_mesh::sweep_management

#endif // CHI_SPDS_H