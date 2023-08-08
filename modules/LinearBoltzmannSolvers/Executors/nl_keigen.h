#ifndef CHITECH_NL_KEIGEN_H
#define CHITECH_NL_KEIGEN_H

#include "physics/SolverBase/chi_solver.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"
#include "A_LBSSolver/IterativeMethods/nl_keigen_ags_solver.h"

#include <petscsnes.h>

namespace lbs
{

class XXNonLinearKEigen : public chi_physics::Solver
{
protected:
  LBSSolver& lbs_solver_;
  std::shared_ptr<NLKEigenAGSContext<Vec,SNES>> nl_context_;
  NLKEigenvalueAGSSolver<Mat,Vec,SNES> nl_solver_;

  bool reinit_phi_1_;
  int num_free_power_its_;

public:
  static chi::InputParameters GetInputParameters();
  explicit XXNonLinearKEigen(const chi::InputParameters& params);

  void Initialize() override;
  void Execute() override;
};

}

#endif // CHITECH_NL_KEIGEN_H
