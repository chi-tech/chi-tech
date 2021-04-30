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
  int num_precursors = 0;
  double k_eff = 1.0;

  // Additional phi vector
  std::vector<double> phi_prev_local;

  // Precursor vector and unknown manager
  std::vector<double> Nj_new_local;
  chi_math::UnknownManager Nj_unk_man;

  // Iterative methods
  void PowerIteration(LBSGroupset& groupset);
  
  // Iterative operations
  void SetKSource(LBSGroupset& groupset,
                  bool apply_mat_src=true,
                  bool suppress_phi_old=false);
  double ComputeProduction();
  void InitializePrecursors();

  // Execute method
  void InitializeKSolver();
  void ExecuteKSolver();
  
};

}

#endif