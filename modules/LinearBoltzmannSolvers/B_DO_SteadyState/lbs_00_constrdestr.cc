#include "lbs_DO_steady_state.h"

//###################################################################
/**Constructor for LBS*/
lbs::DiscOrdSteadyStateSolver::DiscOrdSteadyStateSolver(const std::string& in_text_name) :
  lbs::LBSSolver(in_text_name)
{}

/**Destructor for LBS*/
lbs::DiscOrdSteadyStateSolver::~DiscOrdSteadyStateSolver()
{
  for (auto& groupset : groupsets_)
  {
    CleanUpWGDSA(groupset);
    CleanUpTGDSA(groupset);

    ResetSweepOrderings(groupset);
  }
}