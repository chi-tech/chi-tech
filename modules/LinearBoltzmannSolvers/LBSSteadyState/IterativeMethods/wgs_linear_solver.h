#ifndef CHITECH_WGSLINSOLVEBASE_H
#define CHITECH_WGSLINSOLVEBASE_H

#include "LinearBoltzmannSolvers/LBSSteadyState/Tools/wgs_context.h"

#include "ChiMath/LinearSolver/linear_solver.h"
#include "LinearBoltzmannSolvers/LBSSteadyState/lbs_linear_boltzmann_solver.h"
#include "LinearBoltzmannSolvers/LBSSteadyState/Tools/wgs_context.h"

#include <memory>
#include <vector>
#include <functional>

namespace lbs
{

//################################################################### Class def
/**Linear Solver specialization for Within GroupSet (WGS) solves.*/
template<class MatType, class VecType, class SolverType>
class WGSLinearSolver : public
                        chi_math::LinearSolver<MatType,VecType,SolverType>
{
protected:
  std::vector<double> saved_q_moments_local_;

public:
  typedef std::shared_ptr<WGSContext<MatType,VecType,SolverType>> GSContextPtr;

  /**Constructor.
   * \param gs_context_ptr Context Pointer to abstract context.*/
  explicit WGSLinearSolver(GSContextPtr gs_context_ptr) :
    chi_math::LinearSolver<MatType,VecType,SolverType>
      (IterativeMethodPETScName(gs_context_ptr->groupset_.iterative_method),
       gs_context_ptr)
  {
    auto& groupset = gs_context_ptr->groupset_;
    auto& solver_tol_options = this->ToleranceOptions();
    solver_tol_options.residual_absolute_  = groupset.residual_tolerance;
    solver_tol_options.maximum_iterations_ = groupset.max_iterations;
    solver_tol_options.gmres_restart_interval = groupset.gmres_restart_intvl;
  }

protected:
//  virtual void Setup();
  void PreSetupCallback() override;         //Customized via context
//  virtual void SetOptions();
  void SetSolverContext() override;         //Generic
  void SetConvergenceTest() override;       //Generic
//  virtual void SetMonitor();

  virtual void SetSystemSize() override;    //Customized via context
  virtual void SetSystem() override;        //Generic

  void SetPreconditioner() override;        //Customized via context

  void PostSetupCallback() override;        //Customized via context

//  virtual void Solve();
  void PreSolveCallback() override;         //Customized via context
  void SetRHS() override;                   //Generic + with context elements
  void SetInitialGuess() override;          //Generic
  void PostSolveCallback() override;        //Generic + with context elements
};

}//namespace lbs

#endif //CHITECH_WGSLINSOLVEBASE_H
