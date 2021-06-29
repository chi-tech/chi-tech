#ifndef _k_eigen_solver_h
#define _k_eigen_solver_h

#include "LinearBoltzmannSolver/lbs_linear_boltzmann_solver.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <string>

namespace LinearBoltzmann::KEigenvalue
{

/**A k-eigenvalue neutron transport solver.*/
class Solver : public LinearBoltzmann::Solver
{
private:
  size_t source_event_tag;

public:
  double k_eff = 1.0;

  size_t num_precursors;
  size_t max_num_precursors_per_material;

  chi_math::UnknownManager precursor_uk_man;

  std::vector<double> precursor_new_local;

  // IterativeMethods
  void PowerIteration();

  // Iterative operations
  void SetKSource(LBSGroupset& groupset,
                  std::vector<double>& destination_q,
                  SourceFlags source_flags);
  double ComputeProduction();
  void InitializePrecursors();

  // Execute method
  void InitializeKSolver();
  void ExecuteKSolver();
};

}

#endif