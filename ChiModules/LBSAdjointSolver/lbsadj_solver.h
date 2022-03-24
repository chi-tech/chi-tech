#ifndef LBADJOINTSOLVER_H
#define LBADJOINTSOLVER_H

#include <utility>

#include "LinearBoltzmannSolver/lbs_linear_boltzmann_solver.h"
#include "ChiMath/chi_math.h"

namespace lbs_adjoint
{

struct QOIDesignation
{
  const std::string              name;
  const chi_mesh::LogicalVolume& logical_volume;
  const std::string              lua_functional;

  explicit QOIDesignation(std::string in_name,
                          const chi_mesh::LogicalVolume& in_logical_volume,
                          std::string  in_lua_function_name) :
    name(std::move(in_name)),
    logical_volume(in_logical_volume),
    lua_functional(std::move(in_lua_function_name))
  {}
};

//###################################################################
class AdjointSolver : public lbs::SteadySolver
{
protected:
  std::map<int, std::vector<chi_math::SparseMatrix>> matid_to_S_transpose;

  std::vector<std::pair<QOIDesignation,std::vector<size_t>>> QOI_cell_subscriptions;

public:
  explicit AdjointSolver(const std::string& solver_name);

  void SetSource(LBSGroupset& groupset,
                 std::vector<double>&  destination_q,
                 lbs::SourceFlags source_flags) override;

  void Initialize() override;
  void Execute() override;

  //04
  size_t SetQOI(const std::string& qoi_name,
                const chi_mesh::LogicalVolume& logical_volume,
                const std::string& lua_function_name);
  //05a
  void ExportImportanceMap(const std::string& file_name);
};

}//namespace lbs_adjoint

#endif //LBADJOINTSOLVER_H