#ifndef CHITECH_LBS_ADJOINT_SRC_FUNCTION_H
#define CHITECH_LBS_ADJOINT_SRC_FUNCTION_H

#include "LinearBoltzmannSolvers/A_LBSSolver/SourceFunctions/source_function.h"

namespace lbs
{

class AdjointSourceFunction : public SourceFunction
{
public:
  AdjointSourceFunction(const LBSSolver& lbs_solver,
                        std::shared_ptr<SourceContext>& context);

  void AddAdditionalSources(LBSGroupset& groupset,
                            std::vector<double>& destination_q,
                            const std::vector<double>& phi,
                            SourceFlags source_flags) override
  {
    AddVolumetricQOISources(groupset, destination_q, phi, source_flags);
  }

  void AddVolumetricQOISources(LBSGroupset& groupset,
                               std::vector<double>& destination_q,
                               const std::vector<double>& phi,
                               SourceFlags source_flags);
};


}//namespace lbs

#endif //CHITECH_LBS_ADJOINT_SRC_FUNCTION_H
