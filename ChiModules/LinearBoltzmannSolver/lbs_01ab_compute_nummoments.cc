#include "lbs_linear_boltzmann_solver.h"

//###################################################################
/** Computes the number of moments for the given mesher types*/
void LinearBoltzmann::Solver::ComputeNumberOfMoments()
{
  for (size_t gs = 1; gs < group_sets.size(); ++gs)
    if (group_sets[gs].quadrature->GetMomentToHarmonicsIndexMap()
        != group_sets[0].quadrature->GetMomentToHarmonicsIndexMap())
      throw std::logic_error("LinearBoltzmann::Solver::ComputeNumberOfMoments : "
                             "Moment-to-Harmonics mapping differs between "
                             "groupsets, which is not allowed.");

  num_moments =
    (int)group_sets.front().quadrature->GetMomentToHarmonicsIndexMap().size();

  if (num_moments == 0)
    throw std::logic_error("LinearBoltzmann::Solver::ComputeNumberOfMoments : "
                           "unable to infer number of moments from angular "
                           "quadrature.");
}

