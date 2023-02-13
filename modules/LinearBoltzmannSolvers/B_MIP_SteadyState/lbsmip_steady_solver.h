#ifndef CHITECH_LBSMIP_STEADY_SOLVER_H
#define CHITECH_LBSMIP_STEADY_SOLVER_H

#include "A_LBSSolver/lbs_solver.h"

namespace lbs
{
class MIPSteadyStateSolver : public LBSSolver
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
  void InitializeWGSSolvers() override;
};
}//namespace lbs

#endif //CHITECH_LBSMIP_STEADY_SOLVER_H
