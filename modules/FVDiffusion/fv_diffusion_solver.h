#ifndef FV_DIFFUSION_SOLVER_H
#define FV_DIFFUSION_SOLVER_H

#include "ChiPhysics/SolverBase/chi_solver.h"
#include "ChiMath/PETScUtils/petsc_utils.h"

#include "fv_diffusion_bndry.h"
#include "ChiTimer/chi_timer.h"

#include "ChiConsole/chi_console.h"

#include "ChiMesh/chi_mesh.h"

#include <map>

// forward declaration
namespace chi_mesh
{
  class MeshContinuum;
  typedef std::shared_ptr<MeshContinuum> MeshContinuumPtr;
}
namespace chi_math
{
  class SpatialDiscretization;
  typedef std::shared_ptr<SpatialDiscretization> SDMPtr ;
}

namespace fv_diffusion
{
/** FV diffusion solver
 *
*/
  class Solver : public chi_physics::Solver
  {
  public:
    chi_mesh::MeshContinuumPtr grid_ptr_ = nullptr;

    chi_math::SDMPtr sdm_ptr_ = nullptr;

    size_t num_local_dofs_ = 0;
    size_t num_globl_dofs_ = 0;

    Vec            x_ = nullptr;            // approx solution
    Vec            b_ = nullptr;            // RHS
    Mat            A_ = nullptr;            // linear system matrix

    typedef std::pair<fv_diffusion::BoundaryType,std::vector<double>> BoundaryInfo;
    typedef std::map<uint, BoundaryInfo> BoundaryPreferences;
    BoundaryPreferences      boundary_preferences_;
    std::vector<Boundary>    boundaries_;

    explicit Solver(const std::string& in_solver_name);
    ~Solver() override;

    // void Initialize() override;
    void Initialize() override;
    void Execute() override;

    static double CallLua_iXYZFunction(lua_State* L,
                                       const std::string&,
                                       int,
                                       const chi_mesh::Vector3&);

    void UpdateFieldFunctions();
  };

} // namespace fv_diffusion


#endif //FV_DIFFUSION_SOLVER_H

