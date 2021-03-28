#include "k_eigenvalue_solver.h"

using namespace LinearBoltzmann;

//###################################################################
void KEigenvalue::Solver::InitializeKSolver()
{
  // ----- General LBS init
  Initialize();

  // ----- Number of precursors
  // NOTE: This is only applicable for problems where 
  // precursors only live within one material.
  for (int x=0; x<material_xs.size(); x++) {
    if (material_xs[x]->J > 0) {
      num_precursors = material_xs[x]->J;
      break;
    }
  }

  // ----- Init additional phi vector
  phi_prev_local.resize(phi_old_local.size(),0.0);
}