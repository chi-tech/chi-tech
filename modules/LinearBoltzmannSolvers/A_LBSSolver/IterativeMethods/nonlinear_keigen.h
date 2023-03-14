#ifndef CHITECH_LBS_NONLINEAR_KEIGEN_H
#define CHITECH_LBS_NONLINEAR_KEIGEN_H

namespace lbs
{
  class LBSSolver;

int NonLinearKEigen(LBSSolver& lbs_solver,
                    double nonlinear_abs_tolerance,
                    int    nonlinear_max_iterations,
                    double& k_eff);

}//namespace lbs

#endif //CHITECH_LBS_NONLINEAR_KEIGEN_H
