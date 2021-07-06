#include "k_eigenvalue_solver.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include <chi_log.h>

extern ChiLog& chi_log;


using namespace LinearBoltzmann;

//###################################################################
void KEigenvalue::Solver::InitializeKSolver()
{
  Initialize(); // LinearBoltzmann::Solver initialize

  //======================================== Initialize precursors
  // clear precursor properties
  num_precursors = 0;
  max_num_precursors_per_material = 0;

  // accumulate precursor count and set max
  for (auto& xs : material_xs)
  {
    num_precursors += xs->num_precursors;
    if (xs->num_precursors > max_num_precursors_per_material)
      max_num_precursors_per_material = xs->num_precursors;
  }

  // set use precursors option to false if no precursors
  if (num_precursors == 0)
    options.use_precursors = false;

  // initialize vector and unknown manager if using precursors
  if (options.use_precursors)
  {
    precursor_uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N,
                                max_num_precursors_per_material);

    size_t local_precursor_dof_count =
        discretization->GetNumLocalDOFs(precursor_uk_man);
    precursor_new_local.resize(local_precursor_dof_count, 0.0);
  }//if use precursors
}