#include "k_eigenvalue_solver.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include <chi_log.h>

extern ChiLog& chi_log;


using namespace LinearBoltzmann;

//###################################################################
void KEigenvalue::Solver::InitializeKSolver()
{
  //============================== Parent class initialize
  Initialize();

  //============================== Initialize phi vector
  phi_prev_local.resize(phi_old_local.size(), 0.0);

  //============================== Initialize precursors
  if (options.use_precursors)
  {
    num_precursors = 0;
    precursor_map.clear();

    //============================== Loop over materials
    for (auto& xs : material_xs)
    {
      //============================== Setup map for material
      std::vector<size_t> mat_map(xs->num_precursors);

      //============================== Loop over precursors
      for (size_t j = 0; j < xs->num_precursors; ++j)
        mat_map.emplace_back(num_precursors + j);

      //============================== Increment precursors
      num_precursors += xs->num_precursors;

      //============================== Add mapping to structure
      precursor_map.emplace_back(mat_map);
    }


    if (num_precursors > 0)
    {
      //============================== Initialize precursor unknown manager
      auto pwl = std::static_pointer_cast<SpatialDiscretization_PWLD>(discretization);
      precursor_uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, num_precursors);

      //============================== Initialize precursor vector
      int local_Nj_dof_count = pwl->GetNumLocalDOFs(precursor_uk_man);
      precursor_new_local.resize(local_Nj_dof_count, 0.0);
    }
  }
}