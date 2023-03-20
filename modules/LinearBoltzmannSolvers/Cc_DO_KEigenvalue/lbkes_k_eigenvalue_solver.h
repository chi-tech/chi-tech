#ifndef LBKES_K_EIGENVALUE_SOLVER_H
#define LBKES_K_EIGENVALUE_SOLVER_H

#include "B_DO_Solver/lbs_discrete_ordinates_solver.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <string>

namespace lbs
{

/**A k-eigenvalue solver based on the linear boltzmann transport solver.
\author Zachary Hardy.*/
class DiscOrdKEigenvalueSolver : public lbs::LBSDiscreteOrdinatesSolver
{
private:
  /**The current k-eigenvalue estimate.*/
  double k_eff_ = 1.0;

public:
  explicit DiscOrdKEigenvalueSolver(const std::string& in_text_name) :
    lbs::LBSDiscreteOrdinatesSolver(in_text_name)
  {
    basic_options_.AddOption("K_EIGEN_METHOD", std::string("power"));
    basic_options_.AddOption("K_EIGEN_RESET_SOLUTION", bool(true));
    basic_options_.AddOption("PI_MAX_ITS", int64_t(100));
    basic_options_.AddOption("PI_K_TOL", double(1.0e-10));

    //PISA PowerIteration Synthetic Acceleration
    basic_options_.AddOption("PISA_MIP_L_ABS_TOL", double(1.0e-10));
    basic_options_.AddOption("PISA_MIP_L_MAX_ITS", int64_t(100));

    //PISA NonLinear
    basic_options_.AddOption("PISA_NL_ABS_TOL", double(1.0e-10));
    basic_options_.AddOption("PISA_NL_REL_TOL", double(1.0e-10));
    basic_options_.AddOption("PISA_NL_SOL_TOL", double(1.0e-50));
    basic_options_.AddOption("PISA_NL_MAX_ITS", int64_t(50));

    basic_options_.AddOption("PISA_L_ABS_TOL", double(1.0e-10));
    basic_options_.AddOption("PISA_L_REL_TOL", double(1.0e-10));
    basic_options_.AddOption("PISA_L_MAX_ITS", int64_t(50));

    //PISA PowerIteration
    basic_options_.AddOption("PISA_PI_K_TOL", double(1.0e-10));
    basic_options_.AddOption("PISA_PI_MAX_ITS", int64_t(50));

    basic_options_.AddOption("PISA_VERBOSE_LEVEL", int64_t(0));

    //NLK NonLinear KEigen Value
    basic_options_.AddOption("NLK_ABS_TOL", double(1.0e-8));
    basic_options_.AddOption("NLK_REL_TOL", double(1.0e-8));
    basic_options_.AddOption("NLK_SOL_TOL", double(1.0e-50));
    basic_options_.AddOption("NLK_MAX_ITS", int64_t(50));

    basic_options_.AddOption("NLK_L_REL_TOL", double(1.0e-8));
    basic_options_.AddOption("NLK_L_ABS_TOL", double(1.0e-8));
    basic_options_.AddOption("NLK_L_DIV_TOL", double(1.0e6));
    basic_options_.AddOption("NLK_L_MAX_ITS", int64_t(50));
    basic_options_.AddOption("NLK_GMRES_RESTART_INTVL", int64_t(30));
    basic_options_.AddOption("NLK_GMRES_BRKDN_TOL", double(1.0e6));
  }

  DiscOrdKEigenvalueSolver (const DiscOrdKEigenvalueSolver&) = delete;
  DiscOrdKEigenvalueSolver& operator= (const DiscOrdKEigenvalueSolver&) = delete;

  double GetKeff() const {return k_eff_;}
  void SetKeff(double k_eff) {k_eff_ = k_eff;}

public:
  void Execute() override;
};

}

#endif //LBKES_K_EIGENVALUE_SOLVER_H