#ifndef LBADJOINTSOLVER_H
#define LBADJOINTSOLVER_H

#include <utility>

#include "B_LBSSteadyState/lbs_linear_boltzmann_solver.h"
#include "ChiMath/chi_math.h"

#include "LBSSteadyStateAdjoint/ResponseFunction/lbs_adj_response_function.h"
#include "B_LBSSteadyState/Groupset/lbs_groupset.h"

namespace lbs
{

//###################################################################
class SteadyStateAdjointSolver : public lbs::SteadyStateSolver
{
protected:
  std::map<int, std::vector<chi_math::SparseMatrix>> matid_to_S_transpose;

  typedef std::vector<size_t> VecSize_t;
  typedef std::pair<ResponseFunctionDesignation,VecSize_t> RespFuncAndSubs;
  std::vector<RespFuncAndSubs> response_functions;

public:
  std::vector<std::vector<double>> m_moment_buffers;

public:
  SteadyStateAdjointSolver (const SteadyStateAdjointSolver&) = delete;
  SteadyStateAdjointSolver& operator= (const SteadyStateAdjointSolver&) = delete;

  explicit SteadyStateAdjointSolver(const std::string& solver_name);

  void SetAdjointSource(lbs::LBSGroupset& groupset,
                        std::vector<double>& destination_q,
                        const std::vector<double>& phi,
                        lbs::SourceFlags source_flags);

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

}//namespace lbs

#endif //LBADJOINTSOLVER_H