#ifndef CHITECH_LBS_ADJOINT_SRC_CONTEXT_H
#define CHITECH_LBS_ADJOINT_SRC_CONTEXT_H

#include "LinearBoltzmannSolvers/A_LBSSolver/SourceFunctions/ss_src_context.h"

namespace lbs
{
  class AdjointSourceContext : public lbs::SteadyStateSourceContext
  {
  public:
    AdjointSourceContext() = default;

    double AddSourceMoments() const override;
  };
}//namespace lbs

#endif //CHITECH_LBS_ADJOINT_SRC_CONTEXT_H
