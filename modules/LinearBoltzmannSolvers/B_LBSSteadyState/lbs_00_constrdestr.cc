#include "lbs_linear_boltzmann_solver.h"

//###################################################################
/**Constructor for LBS*/
lbs::SteadyStateSolver::SteadyStateSolver(const std::string& in_text_name) :
  lbs::LBSSolver(in_text_name)
{}

/**Destructor for LBS*/
lbs::SteadyStateSolver::~SteadyStateSolver()
{
  for (auto& groupset : groupsets_)
  {
    CleanUpWGDSA(groupset);
    CleanUpTGDSA(groupset);

    ResetSweepOrderings(groupset);
  }
}