#ifndef CHITECH_LBS_STEADY_STATE_SRC_CONTEXT_H
#define CHITECH_LBS_STEADY_STATE_SRC_CONTEXT_H

#include "source_context.h"

namespace chi_math
{
  class SparseMatrix;
}

namespace lbs
{
  class LBSSolver;
  class LBSGroupset;


class SteadyStateSourceContext : public SourceContext
{
public:
  SteadyStateSourceContext() = default;

  double AddSourceMoments() const override;

  double AddScattering(const chi_math::SparseMatrix& S,
                       const double* phi) const override;

  double AddPromptFission(const std::vector<double>& F_g,
                          const double* phi) const override;

  double AddDelayedFission(const PrecursorList& precursors,
                           const std::vector<double>& nu_delayed_sigma_f,
                           const double* phi) const override;
};

}//namespace lbs

#endif //CHITECH_LBS_STEADY_STATE_SRC_CONTEXT_H
