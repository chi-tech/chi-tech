#ifndef CHI_DIFFUSION_SOLVER_H
#define CHI_DIFFUSION_SOLVER_H

#include "ChiMesh/Cell/cell.h"

#define PROPERTY_D_MAP 5
#define PROPERTY_Q_MAP 6

#include "Modules/DiffusionSolver/Boundaries/chi_diffusion_bndry.h"

#include "Modules/DiffusionSolver/chi_diffusion.h"
#include "ChiPhysics/SolverBase/chi_solver.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

#include "ChiMesh/VolumeMesher/chi_volumemesher.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"
#include "ChiTimer/chi_timer.h"

#include <petscksp.h>

#define PWLC          3
#define PWLD_MIP      4
#define PWLC_GRPS     5
#define PWLD_MIP_GRPS 6
#define PWLD_MIP_GAGG 7

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
  ChiTimer t_assembly;
  ChiTimer t_solve;

  double time_assembly=0.0, time_solve=0.0;
  bool verbose_info=true;
public:
  std::string                              solver_name="Diffusion Solver";
  std::vector<chi_diffusion::Boundary*>    boundaries;
  chi_mesh::MeshContinuumPtr                 grid = nullptr;

  std::shared_ptr<SpatialDiscretization>   discretization;

  chi_math::UnknownManager                 unknown_manager;
  int                                      fem_method = 0;

  int   property_map_D     = 0;
  int   property_map_q     = 1;
  int   property_map_sigma = 2;
  int   material_mode      = DIFFUSION_MATERIALS_REGULAR;
  chi_physics::FieldFunction* D_field     = nullptr;
  chi_physics::FieldFunction* q_field     = nullptr;
  chi_physics::FieldFunction* sigma_field = nullptr;

  bool common_items_initialized=false;

  Vec            x;            // approx solution
  Vec            b;            // RHS
  Mat            A;            // linear system matrix
  KSP            ksp;          // linear solver context
  PC             pc;           // preconditioner context

  PetscReal      norm;         /* norm of solution error */
  PetscErrorCode ierr;         // General error code

  int                            local_dof_count = 0;
  int                            global_dof_count = 0;

  std::vector<double>            pwld_phi_local;

  int    max_iters = 500;
  double residual_tolerance = 1.0e-8;
  int    gi = 0;
  int    G = 1;
  std::string options_string;

public:
  //00
  Solver();
  explicit Solver(std::string in_solver_name);
  virtual ~Solver();
  //01 General
  void GetMaterialProperties(int mat_id,
                             chi_mesh::Cell* cell,
                             int cell_dofs,
                             std::vector<double>& diffCoeff,
                             std::vector<double>& sourceQ,
                             std::vector<double>& sigmaa,
                             int group=0,
                             int moment=0);
  //01a
  void InitializeCommonItems();

  //01b
  int Initialize(bool verbose=true);

  //02a
  int ExecuteS(bool suppress_assembly = false, bool suppress_solve = false);

  void CFEM_Assemble_A_and_b(chi_mesh::Cell& cell, int group=0);

  //02c_c
  void PWLD_Assemble_A_and_b(chi_mesh::Cell& cell,
                             int component=0);
  void PWLD_Assemble_b(chi_mesh::Cell& cell,
                       int component=0);

  //02e_c
  void PWLD_Assemble_A_and_b_GAGG(chi_mesh::Cell& cell);
  void PWLD_Assemble_b_GAGG(chi_mesh::Cell& cell);


  //03b
  double HPerpendicular(chi_mesh::Cell* cell, CellPWLFEValues* fe_view, int f);
  double HPerpendicular(const chi_mesh::Cell& cell,
                        const chi_math::finite_element::UnitIntegralData& fe_intgrl_values,
                        int f);

  static
  uint64_t MapCellLocalNodeIDFromGlobalID(chi_mesh::Cell* cell,
                                          uint64_t node_global_id);
  static
  unsigned int MapCellFace(chi_mesh::Cell* cur_cell,
                           chi_mesh::Cell* adj_cell,
                           unsigned int f);
};

#endif
