#include "lbts_transient_solver.h"

//###################################################################
/**Advances time values.*/
void lbs::DiscOrdTransientSolver::AdvanceTimeValues()
{
  time += dt;
  phi_prev_local = phi_new_local_;
  psi_prev_local = psi_new_local_;
  if (options_.use_precursors)
    precursor_prev_local = precursor_new_local_;
}
