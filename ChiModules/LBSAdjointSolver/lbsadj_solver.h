#ifndef LBADJOINTSOLVER_H
#define LBADJOINTSOLVER_H

#include <utility>

#include "LinearBoltzmannSolver/lbs_linear_boltzmann_solver.h"
#include "ChiMath/chi_math.h"

#include "ResponseFunction/lbs_adj_response_function.h"
#include "LinearBoltzmannSolver/Groupset/lbs_groupset.h"

namespace lbs_adjoint
{

//###################################################################
class AdjointSolver : public lbs::SteadySolver
{
protected:
  std::map<int, std::vector<chi_math::SparseMatrix>> matid_to_S_transpose;

  typedef std::vector<size_t> VecSize_t;
  typedef std::pair<ResponseFunctionDesignation,VecSize_t> RespFuncAndSubs;
  std::vector<RespFuncAndSubs> response_functions;

public:
  AdjointSolver (const AdjointSolver&) = delete;
  AdjointSolver& operator= (const AdjointSolver&) = delete;

  explicit AdjointSolver(const std::string& solver_name);

  void SetSource(lbs::LBSGroupset& groupset,
                 std::vector<double>&  destination_q,
                 lbs::SourceFlags source_flags) override;

  double ComputeInnerProduct();

  void Initialize() override;
  void Execute() override;

  //04
  size_t AddResponseFunction(const std::string& qoi_name,
                             std::shared_ptr<chi_mesh::LogicalVolume> logical_volume,
                             const std::string& lua_function_name);
  //05a
  void ExportImportanceMap(const std::string& file_name);
};

}//namespace lbs_adjoint

#endif //LBADJOINTSOLVER_H