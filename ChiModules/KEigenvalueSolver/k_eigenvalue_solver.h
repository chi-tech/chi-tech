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

  size_t max_iterations = 1000;
  double tolerance = 1.0e-8;

  std::vector<double> phi_prev_local;

  explicit Solver(const std::string& in_text_name) :
    LinearBoltzmann::Solver(in_text_name) {}

  // IterativeMethods
  void PowerIteration();

  // Iterative operations
  void SetKSource(LBSGroupset& groupset,
                  std::vector<double>& destination_q,
                  SourceFlags source_flags);
  double ComputeProduction();
  void InitializePrecursors();

  // Execute method
  void Initialize() override {InitializeKSolver();}
  void InitializeKSolver();
  void ExecuteKSolver();
};

}

#endif