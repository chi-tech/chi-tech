#ifndef DFEM_DIFFUSION_SOLVER_H
#define DFEM_DIFFUSION_SOLVER_H

#include "ChiPhysics/SolverBase/chi_solver.h"
#include "ChiMath/PETScUtils/petsc_utils.h"

#include "dfem_diffusion_bndry.h"
#include "ChiTimer/chi_timer.h"

#include "ChiConsole/chi_console.h"

// forward declaration
namespace chi_mesh
{
class MeshContinuum; 
typedef std::shared_ptr<MeshContinuum> MeshContinuumPtr;
};
namespace chi_math
{
class SpatialDiscretization; 
typedef std::shared_ptr<SpatialDiscretization> SDMPtr ;
};

namespace dfem_diffusion
{
/** DFEM diffusion solver
 * 
*/
class Solver : public chi_physics::Solver
{
private:
  chi_objects::ChiTimer t_assembly;
  chi_objects::ChiTimer t_solve;

  double time_assembly=0.0, time_solve = 0.0;
  bool verbose_info=true;

public:
  chi_mesh::MeshContinuumPtr grid_ptr = nullptr;

  chi_math::SDMPtr sdm_ptr = nullptr;

  size_t num_local_dofs = 0;
  size_t num_globl_dofs = 0;
  chi_math::UnknownManager                 unknown_manager;

  std::vector<double> field;

  Vec            x = nullptr;            // approx solution
  Vec            b = nullptr;            // RHS
  Mat            A = nullptr;            // linear system matrix

  typedef std::pair<BoundaryType,std::vector<double>> BoundaryInfo;
  typedef std::map<uint, BoundaryInfo> BoundaryPreferences;
  BoundaryPreferences      boundary_preferences;
  std::vector<Boundary>   boundaries;

  explicit Solver(const std::string& in_solver_name);
  virtual ~Solver();

  // void Initialize() override;
  void Initialize() override {Initialize(true);}
  void Initialize(bool verbose);

  void Execute() override;

  double HPerpendicular(const chi_mesh::Cell& cell, unsigned int f);

  int MapFaceNodeDisc(const chi_mesh::Cell& cur_cell,
                      const chi_mesh::Cell& adj_cell,
                      const std::vector<chi_mesh::Vector3>& cc_node_locs,
                      const std::vector<chi_mesh::Vector3>& ac_node_locs,
                      size_t ccf, size_t acf,
                      size_t ccfi,
                      double epsilon=1.0e-12);

  static
  double CallLua_iXYZFunction(lua_State* L,
                              const std::string&,
                              const int,
                              const chi_mesh::Vector3&);

};

}; // namespace dfem_diffusion


#endif //DFEM_DIFFUSION_SOLVER_H

