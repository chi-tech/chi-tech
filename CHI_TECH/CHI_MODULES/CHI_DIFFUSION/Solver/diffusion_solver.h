#ifndef _chi_diffusion_solver_h
#define _chi_diffusion_solver_h

#define SOLVER_SYSTEM_EIGEN 1
#define SOLVER_SYSTEM_PETSC 2
#define PROPERTY_D_MAP 5
#define PROPERTY_Q_MAP 6

#include "../Boundaries/chi_diffusion_bndry.h"

#include "../chi_diffusion.h"
#include <CHI_PHYSICS/CHI_SOLVER/chi_solver.h>
#include <spatial_discretization.h>
#include <PiecewiseLinear/pwl.h>

#include <CHI_MESH/CHI_VOLUMEMESHER/chi_volumemesher.h>

#include <CHI_PHYSICS/CHI_FIELDFUNCTION/chi_fieldfunction.h>
#include <ChiTimer/chi_timer.h>

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

struct DIFFUSION_IP_VIEW
{
  int cell_dof_start;
  //std::vector<std::vector<int>> adj_dof_mapping;

  int MapDof(int i)
  {
    return cell_dof_start + i;
  }
};

struct DIFFUSION_IP_BORDERCELL
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

typedef std::vector<std::vector<DIFFUSION_IP_BORDERCELL*>> IP_BORDERCELL_INFO;
typedef std::vector<std::vector<chi_mesh::Cell*>>          IP_BORDERCELLS;
typedef std::vector<std::vector<CellFEView*>>              IP_BORDERFEVIEWS;
typedef std::vector<std::vector<DIFFUSION_IP_VIEW*>>       IP_BORDERIPVIEWS;

//######################################################### Class def
/**Solver for the general diffusion problem.

 <img src="DiffusionMatProp.png" style="width:500px">
 * */
class chi_diffusion::Solver : public chi_physics::Solver
{
private:
  ChiTimer t_assembly;
  ChiTimer t_solve;

  double time_assembly, time_solve;
  bool verbose;
public:
  std::string                              solver_name;
  std::vector<chi_diffusion::Boundary*>    boundaries;
  chi_mesh::MeshContinuum*                 grid;
  SpatialDiscretization*                      discretization;
  SpatialDiscretization_PWL*                  pwl_discr;
  chi_mesh::VolumeMesher*                  mesher;
  int fem_method;

  int   property_map_D;
  int   property_map_q;
  int   property_map_sigma;
  int   material_mode;
  chi_physics::FieldFunction* D_field;
  chi_physics::FieldFunction* q_field;
  chi_physics::FieldFunction* sigma_field;

  bool common_items_initialized;

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
  PetscErrorCode ierr;
  PetscInt       local_rows_from, local_rows_to;

  std::vector<std::vector<int>*> nodal_connections;
  std::vector<std::vector<int>*> nodal_cell_connections;
  std::vector<int>               nodal_boundary_numbers;
  std::vector<int>               nodal_nnz_in_diag;
  std::vector<int>               nodal_nnz_off_diag;

  int                            pwld_local_dof_count;
  int                            pwld_global_dof_count;
  int                            pwld_local_dof_start;
  std::vector<DIFFUSION_IP_VIEW*> ip_cell_views;
  IP_BORDERCELL_INFO             ip_locI_bordercell_info;
  IP_BORDERCELLS                 ip_locI_bordercells;
  IP_BORDERFEVIEWS               ip_locI_borderfeviews;
  IP_BORDERIPVIEWS               ip_locI_borderipviews;

  std::vector<double>            pwld_phi_local;
  std::vector<int>               pwld_cell_dof_array_address;

  int    max_iters;
  double residual_tolerance;
  int    gi;
  int    G;
  std::string options_string;

public:
  //00
  Solver();
  Solver(std::string in_solver_name);
  //01 General
  void GetMaterialProperties(int mat_id, double& diffCoeff,
                                         double& sourceQ,
                                         double& sigmaa);
  void GetMaterialProperties(int mat_id,
                             int cell_glob_index,
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
  //02b_a
  void CFEM_Ab_Slab(int cell_glob_index, chi_mesh::Cell *cell, int group=0);
  //02b_b
  void CFEM_Ab_Polygon(int cell_glob_index, chi_mesh::Cell *cell, int group=0);
  //02b_c
  void CFEM_Ab_Polyhedron(int cell_glob_index, chi_mesh::Cell *cell, int group=0);




  //02c
  int ExecutePWLD_MIP(bool suppress_assembly = false,
                      bool suppress_solve = false);
  //02c_a
  void PWLD_Ab_Slab(int cell_glob_index, chi_mesh::Cell *cell,
                    DIFFUSION_IP_VIEW* cell_ip_view, int group=0);
  void PWLD_b_Slab(int cell_glob_index, chi_mesh::Cell *cell,
                    DIFFUSION_IP_VIEW* cell_ip_view, int group=0);
  //02c_b
  void PWLD_Ab_Polygon(int cell_glob_index, chi_mesh::Cell *cell,
                       DIFFUSION_IP_VIEW* cell_ip_view, int group=0);
  void PWLD_b_Polygon(int cell_glob_index, chi_mesh::Cell *cell,
                       DIFFUSION_IP_VIEW* cell_ip_view, int group=0);
  //02c_c
  void PWLD_Ab_Polyhedron(int cell_glob_index, chi_mesh::Cell *cell,
                          DIFFUSION_IP_VIEW* cell_ip_view, int group=0);
  void PWLD_b_Polyhedron(int cell_glob_index, chi_mesh::Cell *cell,
                          DIFFUSION_IP_VIEW* cell_ip_view, int group=0);

  //02d
  int ExecutePWLD_MIP_GRPS(bool suppress_assembly = false,
                           bool suppress_solve = false);


  //02e
  int ExecutePWLD_MIP_GAGG(bool suppress_assembly = false,
                           bool suppress_solve = false);
  //02e_a
  void PWLD_Ab_Slab_GAGG(int cell_glob_index, chi_mesh::Cell *cell,
                         DIFFUSION_IP_VIEW* cell_ip_view);
  void PWLD_b_Slab_GAGG(int cell_glob_index, chi_mesh::Cell *cell,
                         DIFFUSION_IP_VIEW* cell_ip_view);
  //02e_b
  void PWLD_Ab_Polygon_GAGG(int cell_glob_index, chi_mesh::Cell *cell,
                            DIFFUSION_IP_VIEW* cell_ip_view);
  void PWLD_b_Polygon_GAGG(int cell_glob_index, chi_mesh::Cell *cell,
                            DIFFUSION_IP_VIEW* cell_ip_view);
  //02e_c
  void PWLD_Ab_Polyhedron_GAGG(int cell_glob_index, chi_mesh::Cell *cell,
                               DIFFUSION_IP_VIEW* cell_ip_view);
  void PWLD_b_Polyhedron_GAGG(int cell_glob_index, chi_mesh::Cell *cell,
                              DIFFUSION_IP_VIEW* cell_ip_view);




  //03a
  void ReorderNodesPWLC();

  //03b
  void ReorderNodesPWLD();
  int  MapBorderCell(int locI, int neighbor, int vglob_i);
  void SpawnBorderCell(int locI, int cell_border_index);
  double HPerpendicularPoly(int Nv, double area, double perimeter);
  double HPerpendicularPolyH(int Nf, int Nv, double volume, double area);
  int MapCellDof(chi_mesh::CellSlab* slab_cell, int ig);
  int MapCellDof(chi_mesh::CellPolygon* poly_cell, int ig);
  int MapCellDof(chi_mesh::CellPolyhedron* polyh_cell, int ig);
  int MapCellFace(chi_mesh::CellPolyhedron* polyh_cell,
                  chi_mesh::CellPolyhedron* adjph_cell, int f);

  DIFFUSION_IP_VIEW* GetBorderIPView(int locI, int cell_glob_index);
  CellFEView* GetBorderFEView(int locI, int cell_glob_index);
  chi_mesh::Cell* GetBorderCell(int locI, int cell_glob_index);
};

#endif