#ifndef MG_DIFFUSION_SOLVER_H
#define MG_DIFFUSION_SOLVER_H

#include "ChiPhysics/SolverBase/chi_solver.h"
#include "ChiMath/PETScUtils/petsc_utils.h"

#include "mg_diffusion_bndry.h"
#include "ChiTimer/chi_timer.h"

#include "ChiConsole/chi_console.h"

#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"
#include "ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h"

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

namespace mg_diffusion {

struct KSPAppContext
{
  PetscBool verbose = PETSC_FALSE;
};

struct TwoGridCollapsedInfo
{
  double collapsed_D = 0.0;
  double collapsed_sig_a = 0.0;
  std::vector<double> spectrum;
};

//struct Multigroup_D_and_sigR
//{
//  std::vector<double> Dg;
//  std::vector<double> sigR;
//};

/** Multi-group diffusion solver
 * 
*/
class Solver : public chi_physics::Solver
{
public:
  chi_mesh::MeshContinuumPtr grid_ptr = nullptr;

  chi_math::SDMPtr sdm_ptr = nullptr;

  uint num_groups;
  uint last_fast_group;
  bool do_two_grid = false;

  size_t num_local_dofs = 0;
  size_t num_globl_dofs = 0;

  std::vector<Mat> A;     // linear system matrix for each group
  std::vector<Vec> bext;  // external source vector for each group
  std::vector<Vec> x;     // solution vector for each group
  std::vector<Vec> x_old; // vector of old fluxes

  Vec thermal_dphi = nullptr; // error vector for thermal fluxes
  Vec b = nullptr; // actual rhs vector for the linear system A[g] x[g] = b

  chi_math::PETScUtils::PETScSolverSetup petsc_solver;
  KSPAppContext my_app_context;

  std::vector< std::vector<double> > VF;

//  typedef std::pair<BoundaryType,std::vector<double>> BoundaryInfo;
  typedef std::pair<BoundaryType,std::array<std::vector<double>, 3>>
                                                  BoundaryInfo;

  typedef std::map<uint, BoundaryInfo> BoundaryPreferences;
  BoundaryPreferences     boundary_preferences;
  std::vector<Boundary>   boundaries;

  explicit Solver(const std::string& in_solver_name);
  virtual ~Solver();

  void Initialize() override;

  void Initialize_Materials(std::set<int> &material_ids);
  void Set_BCs(const std::vector<uint64_t>& globl_unique_bndry_ids);
  void Assemble_A_bext();
  void Compute_TwoGrid_Params();
  void Compute_TwoGrid_VolumeFractions();

  void Execute() override;

  void Assemble_RHS(unsigned int g, int64_t iverbose);
  void Assemble_RHS_TwoGrid(int64_t iverbose);
  void SolveOneGroupProblem(unsigned int g, int64_t iverbose);
  void Update_Flux_With_TwoGrid(int64_t iverbose);

  //04
  void UpdateFieldFunctions();

protected:
  std::map<int,std::shared_ptr<chi_physics::TransportCrossSections>>
  matid_to_xs_map;
  std::map<int,std::shared_ptr<chi_physics::IsotropicMultiGrpSource>>
  matid_to_src_map;

  std::map<int, TwoGridCollapsedInfo> map_mat_id_2_tginfo;
//  std::map<int, Multigroup_D_and_sigR> map_mat_id_2_tgXS;

};

} // namespace mg_diffusion


#endif //MG_DIFFUSION_SOLVER_H

