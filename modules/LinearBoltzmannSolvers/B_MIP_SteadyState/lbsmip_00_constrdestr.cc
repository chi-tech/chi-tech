#include "lbsmip_steady_solver.h"

namespace lbs
{

/**Constructor.*/
MIPSteadyStateSolver::MIPSteadyStateSolver(const std::string &in_text_name) :
  LBSSolver(in_text_name)
{

}

MIPSteadyStateSolver::~MIPSteadyStateSolver()
{
  for (auto& groupset : groupsets_)
    CleanUpTGDSA(groupset);
}

}//namespace lbs
