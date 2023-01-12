#ifndef CHI_DIFFUSION_SOLVER_H
#define CHI_DIFFUSION_SOLVER_H

#include "ChiMesh/Cell/cell.h"

#include "DiffusionSolver/Boundaries/chi_diffusion_bndry.h"

#include "DiffusionSolver/chi_diffusion.h"
#include "ChiPhysics/SolverBase/chi_solver.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

#include "ChiMesh/VolumeMesher/chi_volumemesher.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"
#include "ChiTimer/chi_timer.h"


#include <petscksp.h>

#define DIFFUSION_MATERIALS_REGULAR                       10
#define DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTR          11
#define DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF          12
#define DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF_JPART    13
#define DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF_JFULL    14

//######################################################### Class def
/**Solver for the general diffusion problem.

 <img src="DiffusionMatProp.png" style="width:500px">
 * */
class chi_diffusion::Solver : public chi_physics::Solver
{
private:
  chi_objects::ChiTimer t_assembly;
  chi_objects::ChiTimer t_solve;

  double time_assembly=0.0, time_solve=0.0;
  bool verbose_info=true;

public:
  typedef unsigned int uint;
  typedef std::pair<BoundaryType,std::vector<double>> BoundaryInfo;
  typedef std::map<uint, BoundaryInfo> BoundaryPreferences;

public:
  BoundaryPreferences                      boundary_preferences;
  std::vector<chi_diffusion::Boundary*>    boundaries;
  chi_mesh::MeshContinuumPtr                 grid = nullptr;

  std::shared_ptr<chi_math::SpatialDiscretization>   discretization;

  chi_math::UnknownManager                 unknown_manager;

  int   material_mode      = DIFFUSION_MATERIALS_REGULAR;
  std::shared_ptr<chi_physics::FieldFunction> D_field     = nullptr;
  std::shared_ptr<chi_physics::FieldFunction> q_field     = nullptr;
  std::shared_ptr<chi_physics::FieldFunction> sigma_field = nullptr;

  bool common_items_initialized=false;

  Vec            x = nullptr;            // approx solution
  Vec            b = nullptr;            // RHS
  Mat            A = nullptr;            // linear system matrix
  KSP            ksp = nullptr;          // linear solver context
  PC             pc = nullptr;           // preconditioner context

  PetscReal      norm = 0.0;         /* norm of solution error */
  PetscErrorCode ierr = 0;         // General error code

  size_t         local_dof_count = 0;
  size_t         global_dof_count = 0;

  std::vector<double>            pwld_phi_local;

  int    gi = 0;
  int    G = 1;
  std::string options_string;

public:
  //00
  Solver (const Solver&) = delete;
  Solver& operator= (const Solver&) = delete;

  explicit Solver(const std::string& in_solver_name);
  virtual ~Solver();
  //01 General
  void GetMaterialProperties(const chi_mesh::Cell& cell,
                             int cell_dofs,
                             std::vector<double>& diffCoeff,
                             std::vector<double>& sourceQ,
                             std::vector<double>& sigmaa,
                             int group=0,
                             int moment=0);
  //01a
  void InitializeCommonItems();

  //01b
  void Initialize() override {Initialize(true);}
  int Initialize(bool verbose);

  //02a
  void Execute() override
  { ExecuteS(); }
  int ExecuteS(bool suppress_assembly = false, bool suppress_solve = false);

  void CFEM_Assemble_A_and_b(chi_mesh::Cell& cell, int group=0);

  //02c_c
  void PWLD_Assemble_A_and_b(const chi_mesh::Cell& cell,
                             int component=0);
  void PWLD_Assemble_b(const chi_mesh::Cell& cell,
                       int component=0);

  //02e_c
//  void PWLD_Assemble_A_and_b_GAGG(const chi_mesh::Cell& cell);
//  void PWLD_Assemble_b_GAGG(const chi_mesh::Cell& cell);


  //03b
  double HPerpendicular(const chi_mesh::Cell& cell,
                        const chi_math::finite_element::UnitIntegralData& fe_intgrl_values,
                        unsigned int f);

  static
  uint64_t MapCellLocalNodeIDFromGlobalID(const chi_mesh::Cell& cell,
                                          uint64_t node_global_id);

  static
  unsigned int MapCellFace(const chi_mesh::Cell& cur_cell,
                           const chi_mesh::Cell& adj_cell,
                           unsigned int f);

  void UpdateFieldFunctions();
};

#endif
