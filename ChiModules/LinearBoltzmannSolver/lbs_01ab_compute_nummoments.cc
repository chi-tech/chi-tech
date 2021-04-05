#include "lbs_linear_boltzmann_solver.h"

//###################################################################
/** Computes the number of moments for the given mesher types*/
void LinearBoltzmann::Solver::ComputeNumberOfMoments()
{
  num_moments =
    group_sets.front().quadrature->GetMomentToHarmonicsIndexMap().size();

  if (num_moments == 0)
    throw std::logic_error("LinearBoltzmann::Solver::ComputeNumberOfMoments : "
                           "unable to infer number of moments from angular "
                           "quadrature.");
}

