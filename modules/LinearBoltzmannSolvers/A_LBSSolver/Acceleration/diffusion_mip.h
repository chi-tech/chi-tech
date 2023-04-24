#ifndef CHITECH_LBS_DIFFUSION_MIP_H
#define CHITECH_LBS_DIFFUSION_MIP_H

#include "acceleration.h"
#include "ChiMath/UnknownManager/unknown_manager.h"
#include "petscksp.h"
#include "ChiLua/chi_lua.h"



//############################################### Forward declarations
namespace chi_mesh
{
  class MeshContinuum;
  class Cell;
  struct Vector3;
}

namespace chi_math
{
  class SpatialDiscretization;
}

namespace lbs
{
  struct UnitCellMatrices;
}

//############################################### Namespace lbs::acceleration
namespace lbs::acceleration
{

struct Multigroup_D_and_sigR;



/**Generalized diffusion solver for both WGDSA and TGDSA based on the MIP-method
 * of Bruno Turcksin and Jean Ragusa.*/
class DiffusionMIPSolver
{
protected:
  typedef std::map<int,Multigroup_D_and_sigR> MatID2XSMap;
protected:
  const std::string text_name_;
  const chi_mesh::MeshContinuum& grid_;
  const chi_math::SpatialDiscretization& sdm_;
  const chi_math::UnknownManager uk_man_;

  const std::map<uint64_t, BoundaryCondition> bcs_;

  const MatID2XSMap mat_id_2_xs_map_;

  const std::vector<UnitCellMatrices>& unit_cell_matrices_;

  const int64_t num_local_dofs_;
  const int64_t num_global_dofs_;

  Mat A_ = nullptr;
  Vec rhs_ = nullptr;
  KSP ksp_ = nullptr;

public:
  struct Options
  {
    double      residual_tolerance = 1.0e-4; ///< Residual tol. relative to rhs
    int         max_iters = 100;             ///< Maximum iterations
    bool        verbose = false;             ///< Verbosity flag
    bool        perform_symmetry_check = false; ///< For debugging only (very expensive)
    std::string source_lua_function;       ///<for mms
    std::string ref_solution_lua_function; ///<for mms
    std::string additional_options_string;
    double      penalty_factor = 4.0;
  }options;
public:
  const chi_math::UnknownManager& UnknownStructure() const
  {
    return uk_man_;
  }
  //00
  DiffusionMIPSolver(std::string text_name,
                     const chi_math::SpatialDiscretization& sdm,
                     const chi_math::UnknownManager& uk_man,
                     std::map<uint64_t, BoundaryCondition> bcs,
                     MatID2XSMap map_mat_id_2_xs,
                     const std::vector<UnitCellMatrices>& unit_cell_matrices,
                     bool verbose=false);
  std::string TextName() const {return text_name_;}
  const Vec& RHS() const;
  //00a
  void Initialize();

  //01a
  void AssembleAand_b_wQpoints(const std::vector<double>& q_vector);
  //01b
  void Assemble_b_wQpoints(const std::vector<double>& q_vector);

  //01c
  void AssembleAand_b(const std::vector<double>& q_vector);
  //01d
  void Assemble_b(const std::vector<double>& q_vector);
  void Assemble_b(Vec petsc_q_vector);

  //02
  void Solve(std::vector<double>& solution, bool use_initial_guess=false);
  void Solve(Vec petsc_solution, bool use_initial_guess=false);
  //05
  double HPerpendicular(const chi_mesh::Cell& cell, unsigned int f);

  int MapFaceNodeDisc(const chi_mesh::Cell& cur_cell,
                      const chi_mesh::Cell& adj_cell,
                      const std::vector<chi_mesh::Vector3>& cc_node_locs,
                      const std::vector<chi_mesh::Vector3>& ac_node_locs,
                      size_t ccf, size_t acf,
                      size_t ccfi,
                      double epsilon=1.0e-12);
  static
  double CallLuaXYZFunction(lua_State* L, const std::string& lua_func_name,
                            const chi_mesh::Vector3& xyz);

  std::pair<size_t, size_t> GetNumPhiIterativeUnknowns();

  virtual ~DiffusionMIPSolver();
};

}//namespace lbs::acceleration

#endif //CHITECH_LBS_DIFFUSION_MIP_H
