#ifndef CHITECH_LBS_STEADY_STATE_H
#define CHITECH_LBS_STEADY_STATE_H

#include "LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

namespace lbs
{

class SteadyStateSolver : public chi_physics::Solver
{
protected:
  LBSSolver& lbs_solver_;

public:
  static chi_objects::InputParameters GetInputParameters();

  explicit SteadyStateSolver(const chi_objects::InputParameters& params);

  void Initialize() override;
  void Execute() override;
};

}

#endif // CHITECH_LBS_STEADY_STATE_H
