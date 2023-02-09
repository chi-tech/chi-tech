#include "lbts_transient_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

/**Constructor.*/
lbs::TransientSolver::TransientSolver(const std::string &in_text_name) :
  lbs::KEigenvalueSolver(in_text_name)
{
  chi::log.Log() << TextName() << " created.";
}

/**Destructor*/
lbs::TransientSolver::~TransientSolver()
{

}