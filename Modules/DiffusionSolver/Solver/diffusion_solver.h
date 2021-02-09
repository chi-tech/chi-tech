#ifndef _chi_diffusion_solver_h
#define _chi_diffusion_solver_h

#include "ChiMesh/Cell/cell.h"

#define SOLVER_SYSTEM_EIGEN 1
#define SOLVER_SYSTEM_PETSC 2
#define PROPERTY_D_MAP 5
#define PROPERTY_Q_MAP 6

#include "Modules/DiffusionSolver/Boundaries/chi_diffusion_bndry.h"

#include "Modules/DiffusionSolver/chi_diffusion.h"
#include "ChiPhysics/SolverBase/chi_solver.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/SpatialDiscretization/PiecewiseLinear/pwl.h"

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

struct DiffusionIPCellView
{
  int cell_dof_start;
  //std::vector<std::vector<int>> adj_dof_mapping;

  int MapDof(int i)
  {
    return cell_dof_start + i;
  }
};

struct DiffusionIPBorderCell
{
  int cell_glob_index;
  int cell_dof_start;
  int cell_type;
  int cell_mat_id;
  int cell_dof_count;
  int cell_face_count;
  std::vector<int> v_indices;

  std::vector<std::vector<int>> face_v_indices;
};

typedef std::vector<std::vector<DiffusionIPBorderCell*>> IP_BORDERCELL_INFO;
typedef std::vector<std::vector<chi_mesh::Cell*>>          IP_BORDERCELLS;
typedef std::vector<std::vector<CellPWLFEValues*>>              IP_BORDERFEVIEWS;
typedef std::vector<std::vector<DiffusionIPCellView*>>       IP_BORDERIPVIEWS;

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

  std::shared_ptr<SpatialDiscretization>     discretization;
  std::shared_ptr<SpatialDiscretization_PWL> pwl_sdm;

  chi_math::NodalVariableStructure                 unknown_manager;
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

  Vec*           xg;           // Groupwise approx solution
  Vec*           bg;           // Groupwise RHS
  Mat*           Ag;           // Groupwise linear system matrix
  KSP*           kspg;         // Groupwise linear solver context
  PC*            pcg;          // Groupwise preconditioner context

  Vec            xref;         // Reference solution vector
  Vec            bref;         // Reference RHS
  Mat            Aref;         // Reference linear system matrix

  PetscReal      norm;         /* norm of solution error */
  PetscErrorCode ierr;         // General error code

  std::vector<std::vector<int>>  nodal_connections;
  std::vector<std::vector<int>>  nodal_cell_connections;
  std::vector<int>               nodal_boundary_numbers;
  std::vector<int>               nodal_nnz_in_diag;
  std::vector<int>               nodal_nnz_off_diag;

  int                            local_dof_count = 0;
  int                            global_dof_count = 0;
  int                            pwld_local_dof_start = 0;
  std::vector<DiffusionIPCellView*> ip_cell_views;
  IP_BORDERCELL_INFO             ip_locI_bordercell_info;
  IP_BORDERCELLS                 ip_locI_bordercells;
  IP_BORDERFEVIEWS               ip_locI_borderfeviews;
  IP_BORDERIPVIEWS               ip_locI_borderipviews;

  std::vector<double>            pwld_phi_local;
  std::vector<int>               pwld_cell_dof_array_address;

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
//  void GetMaterialProperties(int mat_id, double& diffCoeff,
//                                         double& sourceQ,
//                                         double& sigmaa);
  void GetMaterialProperties(int mat_id,
                             chi_mesh::Cell* cell,
                             int cell_dofs,
                             std::vector<double>& diffCoeff,
                             std::vector<double>& sourceQ,
                             std::vector<double>& sigmaa,
                             int group=0,
                             int moment=0);
  bool ApplyDirichletI(int ir, int* ir_boundary_type,
                       int iref=-1,
                       bool suppress_assembly=false);
  bool ApplyDirichletJ(int jr, int ir,
                       double jr_mat_entry,
                       int* ir_boundary_type, int iref=-1);
  //01a
  void InitializeCommonItems();
  void PWLDBuildSparsityPattern();
  void PWLCBuildSparsityPattern();
  //01b
  int Initialize(bool verbose=true);
  //01c
  int InitializePWLC(bool verbose=true);
  //01d
  int InitializePWLD(bool verbose=true);
  //01e
  int InitializePWLDGroups(bool verbose=true);
  //01f
  int InitializePWLDGrpAgg(bool verbose=true);


  //02a
  int ExecuteS(bool suppress_assembly = false, bool suppress_solve = false);



  //02b
  int  ExecutePWLC(bool suppress_assembly = false,
                   bool suppress_solve = false);
  void CFEM_Assemble_A_and_b(int cell_glob_index, chi_mesh::Cell *cell, int group=0);




  //02c
  int ExecutePWLD_MIP(bool suppress_assembly = false,
                      bool suppress_solve = false);
  //02c_c
  void PWLD_Assemble_A_and_b(int cell_glob_index,
                             chi_mesh::Cell *cell,
                             DiffusionIPCellView* cell_ip_view,
                             int component=0,
                             int component_block_offset=1);
  void PWLD_Assemble_b(int cell_glob_index,
                       chi_mesh::Cell *cell,
                       DiffusionIPCellView* cell_ip_view,
                       int component=0,
                       int component_block_offset=1);


  //02d
  int ExecutePWLD_MIP_GRPS(bool suppress_assembly = false,
                           bool suppress_solve = false);

  //02e
  int ExecutePWLD_MIP_GAGG(bool suppress_assembly = false,
                           bool suppress_solve = false);

  //02e_c
  void PWLD_Assemble_A_and_b_GAGG(int cell_glob_index, chi_mesh::Cell *cell,
                               DiffusionIPCellView* cell_ip_view);
  void PWLD_Assemble_b_GAGG(int cell_glob_index, chi_mesh::Cell *cell,
                              DiffusionIPCellView* cell_ip_view);




  //03a
  void ReorderNodesPWLC();

  //03b
  void ReorderNodesPWLD();
  int  MapBorderCell(int locI, int neighbor, int vglob_i);
  void SpawnBorderCell(int locI, int cell_border_index);


  double HPerpendicular(chi_mesh::Cell* cell, CellPWLFEValues* fe_view, int f);
  int MapCellDof(chi_mesh::Cell* cell, int ig);
  int MapCellFace(chi_mesh::Cell* cur_cell,
                  chi_mesh::Cell* adj_cell, int f);

  DiffusionIPCellView* GetBorderIPView(int locI, int cell_glob_index);
  CellPWLFEValues* GetBorderFEView(int locI, int cell_glob_index);
  chi_mesh::Cell* GetBorderCell(int locI, int cell_glob_index);
};

#endif
