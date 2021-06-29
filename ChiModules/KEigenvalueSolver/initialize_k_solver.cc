#include "k_eigenvalue_solver.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include <chi_log.h>

extern ChiLog& chi_log;


using namespace LinearBoltzmann;

//###################################################################
void KEigenvalue::Solver::InitializeKSolver()
{
  Initialize(); // LinearBoltzmann::Solver initialize

  phi_prev_local.resize(phi_old_local.size(), 0.0);

  //======================================== Initialize precursors
  if (options.use_precursors)
  {
    //clear precursor properties
    num_precursors = 0.0;
    max_num_precursors_per_material = 0.0;

    //accumulate precursor count and set max
    for (auto& xs : material_xs)
    {
      num_precursors += xs->num_precursors;
      if (xs->num_precursors > max_num_precursors_per_material)
        max_num_precursors_per_material = xs->num_precursors;
    }


    if (num_precursors > 0)
    {
      typedef SpatialDiscretization_PWLD  PWLD;
      auto pwl = std::static_pointer_cast<PWLD>(discretization);

      precursor_uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N,
                                  max_num_precursors_per_material);

      int local_precursor_dof_count = pwl->GetNumLocalDOFs(precursor_uk_man);
      precursor_new_local.resize(local_precursor_dof_count, 0.0);
    }
  }//if use precursors
}