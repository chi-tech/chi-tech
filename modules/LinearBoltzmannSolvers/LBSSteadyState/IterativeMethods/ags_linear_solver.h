#ifndef CHITECH_AGS_LINEAR_SOLVER_H
#define CHITECH_AGS_LINEAR_SOLVER_H

#include "ChiMath/LinearSolver/linear_solver.h"
#include "LBSSteadyState/Tools/ags_context.h"

namespace lbs
{

//################################################################### Class def
/**Linear Solver specialization for Within GroupSet (WGS) solves.*/
template<class MatType, class VecType, class SolverType>
class AGSLinearSolver : public
                        chi_math::LinearSolver<MatType,VecType,SolverType>
{
protected:
  int groupspan_first_id_ = 0;
  int groupspan_last_id_ = 0;
public:
  typedef std::shared_ptr<AGSContext<MatType,VecType,SolverType>> AGSContextPtr;

  /**Constructor.
   * \param gs_context_ptr Context Pointer to abstract context.*/
  explicit AGSLinearSolver(std::string iterative_method,
                           AGSContextPtr ags_context_ptr,
                           int groupspan_first_id,
                           int groupspan_last_id) :
    chi_math::LinearSolver<MatType,VecType,SolverType>
      (std::move(iterative_method),ags_context_ptr),
    groupspan_first_id_(groupspan_first_id),
    groupspan_last_id_(groupspan_last_id)
  {}

  int GroupSpanFirstID() const {return groupspan_first_id_;}
  int GroupSpanLastID() const {return groupspan_last_id_;}

protected:
  /*void PreSetupCallback() override;  */   //Customized via context
  /*virtual void SetOptions();         */
  /*void SetSolverContext() override;  */   //Generic
  /*void SetConvergenceTest() override;*/   //Generic
  /*virtual void SetMonitor();         */

  virtual void SetSystemSize() override;    //Customized via context
  virtual void SetSystem() override;        //Generic
  void SetPreconditioner() override;        //Customized via context

  /*void PostSetupCallback() override;*/    //Customized via context

public:
  /*virtual void Setup();*/

protected:
  /*void PreSolveCallback() override;*/     //Customized via context
  void SetRHS() override;                   //Generic + with context elements
  void SetInitialGuess() override;          //Generic
  /*void PostSolveCallback() override;*/    //Generic + with context elements
public:
  void Solve() override;

public:

};

}//namespace lbs

#endif //CHITECH_AGS_LINEAR_SOLVER_H
