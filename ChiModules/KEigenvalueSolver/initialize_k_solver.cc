#include "k_eigenvalue_solver.h"

using namespace LinearBoltzmann;

//###################################################################
void KEigenvalue::Solver::InitializeKSolver()
{
  // ----- General LBS init
  Initialize();

  // ----- Number of precursors
  /** NOTE: This computes the total number of precursors
   //       in the problem across all materials and assigns
   //       a global mapping to the material. For example,
   //       if material 0 has 6 precursors and material 1
   //       has 6, the global mapping will give material 0
   //       precursor IDs 0-5 and material 1 will receive
   //       precursor IDs 6-11.
  **/
  num_precursors = 0;
  for (int x = 0; x < material_xs.size(); x++) {
    if (material_xs[x]->J > 0) {
      // Define the precursor mapping for this material
      for (int j = 0; j < material_xs[x]->J; ++j) {
        material_xs[x]->precursor_map[j] = num_precursors + j;
      }

      // Increment the total number of precursors
      num_precursors += material_xs[x]->J;
    }
  }

  // ----- Init additional phi vector
  phi_prev_local.resize(phi_old_local.size(),0.0);
}