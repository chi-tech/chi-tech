#ifndef _k_eigen_solver_h
#define _k_eigen_solver_h

#include "LinearBoltzmannSolver/lbs_linear_boltzmann_solver.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <string>

namespace  LinearBoltzmann::KEigenvalue
{

/**A k-eigenvalue neutron transport solver.*/
class Solver : public LinearBoltzmann::Solver
{
private:
  size_t source_event_tag;

public:
  double k_eff = 1.0;

  size_t num_precursors;

  std::vector<double> phi_prev_local;

  // Precursor vector and unknown manager
  std::vector<double> precursor_new_local;
  chi_math::UnknownManager precursor_uk_man;

  // This structure maps local precursor indices to
  // global precursor indices. This is used to assign
  // a global numbering system for precursors which
  // becomes important when multiple materials with
  // precursors exist within a problem.
  std::vector<std::vector<size_t>> precursor_map;

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