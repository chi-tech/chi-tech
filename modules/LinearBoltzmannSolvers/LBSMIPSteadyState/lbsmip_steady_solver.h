#ifndef CHITECH_LBSMIP_STEADY_SOLVER_H
#define CHITECH_LBSMIP_STEADY_SOLVER_H

#include "LBSSteadyState/lbs_linear_boltzmann_solver.h"

namespace lbs
{
class MIPSteadyStateSolver : public SteadyStateSolver
{
protected:
  typedef std::shared_ptr<acceleration::DiffusionMIPSolver> MIPSolverPtr;
public:
  std::vector<MIPSolverPtr> gs_mip_solvers_;
public:
  //00
  explicit MIPSteadyStateSolver(const std::string& in_text_name);
  //01
  void Initialize() override;
  //02
  void Execute() override;
  void SolveGroupset(LBSGroupset& groupset) override;
};
}//namespace lbs

#endif //CHITECH_LBSMIP_STEADY_SOLVER_H
