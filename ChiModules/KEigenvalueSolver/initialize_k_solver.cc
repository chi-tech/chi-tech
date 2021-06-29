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
  if (options.use_precursors)
  {
    //clear precursor properties
    num_precursors = 0.0;
    max_num_precursors_per_node = 0.0;

    //========== Loop over materials
    for (auto& xs : material_xs)
    {
      num_precursors += xs->num_precursors;
      if (xs->num_precursors > max_num_precursors_per_node)
        max_num_precursors_per_node = xs->num_precursors;
    }

    // ----- Initialize precursor unknown manager and vector
    if (num_precursors > 0)
    {
      typedef SpatialDiscretization_PWLD  PWLD;
      auto pwl = std::static_pointer_cast<PWLD>(discretization);

      Nj_unk_man.AddUnknown(chi_math::UnknownType::VECTOR_N,
                            max_num_precursors_per_node);

      int local_Nj_dof_count = pwl->GetNumLocalDOFs(Nj_unk_man);
      Nj_new_local.resize(local_Nj_dof_count, 0.0);
    }
  }
}