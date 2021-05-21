#include "k_eigenvalue_solver.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include <chi_log.h>
extern ChiLog& chi_log;


using namespace LinearBoltzmann;

//###################################################################
void KEigenvalue::Solver::InitializeKSolver()
{
  // ----- General LBS init
  Initialize();

  // ----- Init additional phi vector
  phi_prev_local.resize(phi_old_local.size(),0.0);

  // ----- Init precursors
  if (options.use_precursors) {
    // ----- Count precursors and define mapping
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
      if (material_xs[x]->num_precursors > 0) {
        // Define the precursor mapping for this material
        for (int j = 0; j < material_xs[x]->num_precursors; ++j) {
          material_xs[x]->precursor_map[j] = num_precursors + j;
        }

        // Increment the total number of precursors
        num_precursors += material_xs[x]->num_precursors;
      }
    }

    // -----Initialize precursor unknown manager and vector
    if (num_precursors > 0) {
      auto pwl = std::static_pointer_cast<SpatialDiscretization_PWLD>(discretization);
      Nj_unk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, num_precursors);

      int local_Nj_dof_count = pwl->GetNumLocalDOFs(Nj_unk_man);
      Nj_new_local.resize(local_Nj_dof_count, 0.0);
    }
  }
}