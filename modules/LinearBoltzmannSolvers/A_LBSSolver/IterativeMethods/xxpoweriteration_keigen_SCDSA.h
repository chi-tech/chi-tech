#ifndef CHITECH_XXPOWERITERATION_KEIGEN_SCDSA_H
#define CHITECH_XXPOWERITERATION_KEIGEN_SCDSA_H

#include "xxpoweriteration_keigen.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Acceleration/diffusion.h"

namespace lbs
{

class XXPowerIterationKEigenSCDSA : public XXPowerIterationKEigen
{
  typedef std::shared_ptr<acceleration::DiffusionSolver> DiffusionSolverPtr;

protected:
  int accel_pi_max_its_;
  double accel_pi_k_tol_;
  bool accel_pi_verbose_;
  DiffusionSolverPtr diffusion_solver_ = nullptr;

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
};

} // namespace lbs

#endif // CHITECH_XXPOWERITERATION_KEIGEN_SCDSA_H
