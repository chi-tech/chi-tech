#ifndef CHITECH_LBS_TRANSIENT_H
#define CHITECH_LBS_TRANSIENT_H

#include "LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"
namespace chi_math
{
class TimeIntegration;
}

namespace lbs
{

class TransientSolver : public chi_physics::Solver
{
protected:
  LBSSolver& lbs_solver_;
  std::shared_ptr<chi_math::TimeIntegration> time_integration_;

public:
  static chi_objects::InputParameters GetInputParameters();
  explicit TransientSolver(const chi_objects::InputParameters& params);

  void Initialize() override;
  void Execute() override;
  void Step() override;
  void Advance() override;
};

}

#endif // CHITECH_LBS_TRANSIENT_H
