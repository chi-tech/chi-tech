#ifndef LINEAR_BOLTZMANN_SOLVER_H
#define LINEAR_BOLTZMANN_SOLVER_H

#include "B_DO_Solver/lbs_discrete_ordinates_solver.h"

namespace lbs
{
//################################################################### Class def
/**A neutral particle transport solver.*/
class DiscOrdSteadyStateSolver : public DiscreteOrdinatesSolver
{
public:
  static chi_objects::InputParameters GetInputParameters();

  explicit DiscOrdSteadyStateSolver(const chi_objects::InputParameters& params);

  explicit DiscOrdSteadyStateSolver(const std::string& in_text_name) :
    lbs::DiscreteOrdinatesSolver(in_text_name) {}

  DiscOrdSteadyStateSolver (const DiscOrdSteadyStateSolver&) = delete;
  DiscOrdSteadyStateSolver& operator= (const DiscOrdSteadyStateSolver&) = delete;

public:
  //02
  void Execute() override;

};

}

#endif
