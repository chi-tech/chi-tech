#ifndef CHITECH_LBTS_TRANSIENT_SOLVER_H
#define CHITECH_LBTS_TRANSIENT_SOLVER_H

#include "LBKEigenvalueSolver/lbkes_k_eigenvalue_solver.h"
#include "ChiMath/chi_math_time_stepping.h"

typedef chi_mesh::sweep_management::SweepChunk SweepChunk;

namespace lbs
{

//################################################################### Class def
/**A transient neutral particle transport solver.
 *
\author Zachary Hardy.*/
class TransientSolver : public KEigenvalueSolver
{
public:
  chi_math::SteppingMethod method = chi_math::SteppingMethod::CRANK_NICHOLSON;
  struct Options
  {
    int verbosity_level = 1;
    bool inhibit_advance = false;
    double t_final = 0.1;
    int max_time_steps = 10;
    std::string console_call_back_function;
  }transient_options;

  /**Temporal domain and discretization information.*/
  double dt = 2.0e-3;
  double time = 0.0;

protected:
  /**Previous time step vectors.*/
  std::vector<double> phi_prev_local;
  std::vector<double> precursor_prev_local;
  std::vector<std::vector<double>> psi_prev_local;

  /**Fission rate vector*/
  std::vector<double> fission_rate_local;

public:
  explicit TransientSolver(const std::string& in_text_name);

  //01
  void Initialize() override;
  //02
  void Execute() override;
  void Step() override;
  void AdvanceTimeValues();

  //Iterative operations
  std::shared_ptr<SweepChunk> SetTransientSweepChunk(LBSGroupset& groupset);
  void SetTransientSource(LBSGroupset& groupset,
                          std::vector<double>&  destination_q,
                          SourceFlags source_flags);
  double ComputeFissionRate(bool previous) override;
  double ComputeBeta();
  void   PostStepCallBackFunction() const;

  //precursors
  void StepPrecursors();

  virtual ~TransientSolver();
};

}//namespace lbs

#endif //CHITECH_LBTS_TRANSIENT_SOLVER_H
