#ifndef CHITECH_PI_KEIGEN_SCDSA_H
#define CHITECH_PI_KEIGEN_SCDSA_H

#include "pi_keigen.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Acceleration/diffusion.h"

namespace chi_math
{
class VectorGhostCommunicator;
}

namespace lbs
{

class XXPowerIterationKEigenSCDSA : public XXPowerIterationKEigen
{
  typedef std::shared_ptr<acceleration::DiffusionSolver> DiffusionSolverPtr;
  typedef std::shared_ptr<chi_math::VectorGhostCommunicator> VecGhostCommPtr;
  typedef std::shared_ptr<chi_math::SpatialDiscretization> SDMPtr;

protected:
  int accel_pi_max_its_;
  double accel_pi_k_tol_;
  bool accel_pi_verbose_;
  DiffusionSolverPtr diffusion_solver_ = nullptr;

  const std::string diffusion_solver_sdm_;

  SDMPtr continuous_sdm_ptr_ = nullptr;
  bool requires_ghosts_ = false;
  struct GhostInfo
  {
    VecGhostCommPtr vector_ghost_communicator = nullptr;
    std::map<int64_t, int64_t> ghost_global_id_2_local_map;
  };
  GhostInfo pwld_ghost_info_;

public:
  static chi_objects::InputParameters GetInputParameters();
  explicit XXPowerIterationKEigenSCDSA(
    const chi_objects::InputParameters& params);

  void Execute() override;

  std::vector<double> CopyOnlyPhi0(const LBSGroupset& groupset,
                                   const std::vector<double>& phi_in);

  void ProjectBackPhi0(const LBSGroupset& groupset,
                       const std::vector<double>& input,
                       std::vector<double>& output);

  GhostInfo MakePWLDVecGhostCommInfo();
  std::vector<double> MakeContinuousVersion(const std::vector<double>& input);
};

} // namespace lbs

#endif // CHITECH_PI_KEIGEN_SCDSA_H