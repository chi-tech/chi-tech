#ifndef CHITECH_LBS_TD_SRC_CONTEXT_H
#define CHITECH_LBS_TD_SRC_CONTEXT_H

#include "LinearBoltzmannSolvers/A_LBSSolver/SourceFunctions/ss_src_context.h"

#include "ChiMath/chi_math_time_stepping.h"

namespace lbs
{
  class TransientSourceContext : SteadyStateSourceContext
  {
  private:
    double& dt_;
    chi_math::SteppingMethod& method_;
  public:
    TransientSourceContext(double& ref_dt,
                           chi_math::SteppingMethod& method) :
                           dt_(ref_dt),
                           method_(method)
    {}

    double AddDelayedFission(const PrecursorList& precursors,
                             const std::vector<double>& nu_delayed_sigma_f,
                             const double* phi) const override;
  };
}

#endif //CHITECH_LBS_TD_SRC_CONTEXT_H
