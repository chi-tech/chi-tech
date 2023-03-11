#ifndef CHITECH_LBS_SOURCE_FUNCTION_H
#define CHITECH_LBS_SOURCE_FUNCTION_H

#include "LinearBoltzmannSolvers/A_LBSSolver/lbs_structs.h"

#include <memory>
#include <utility>

namespace chi_mesh
{
  class MeshContinuum;
  class Cell;
}
typedef std::shared_ptr<const chi_mesh::MeshContinuum> ConstGridPtr;

namespace chi_math
{
  class AngularQuadrature;
}

namespace lbs
{
class LBSSolver;
class LBSGroupset;

class SourceContext;

//###################################################################
/**Implements a customizable source function using virtual methods and
 * an additionally customizable source context.*/
class SourceFunction
{
protected:
  const LBSSolver& lbs_solver_;
  const chi_mesh::MeshContinuum& grid_;
  std::vector<chi_math::AngularQuadrature::HarmonicIndices> m_to_ell_em_map_;
  std::shared_ptr<SourceContext> context_;

public:
  SourceFunction(const LBSSolver& lbs_solver,
                 std::shared_ptr<SourceContext>& context);

  virtual void operator()(LBSGroupset& groupset,
                          std::vector<double>& destination_q,
                          const std::vector<double>& phi,
                          SourceFlags source_flags);

  virtual void AddAdditionalSources(LBSGroupset& groupset,
                                    std::vector<double>& destination_q,
                                    const std::vector<double>& phi,
                                    SourceFlags source_flags)
  {
    AddPointSources(groupset, destination_q, phi, source_flags);
  }

  void AddPointSources(LBSGroupset& groupset,
                       std::vector<double>& destination_q,
                       const std::vector<double>& phi,
                       SourceFlags source_flags);
};

}//namespace lbs

#endif //CHITECH_LBS_SOURCE_FUNCTION_H
