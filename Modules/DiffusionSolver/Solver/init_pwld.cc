#include "diffusion_solver.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiTimer/chi_timer.h>

#include <chi_mpi.h>
#include <chi_log.h>
#include <ChiPhysics/chi_physics.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;
extern ChiPhysics chi_physics_handler;

PetscErrorCode
DiffusionConvergenceTestNPT(KSP ksp, PetscInt n, PetscReal rnorm,
                            KSPConvergedReason* convergedReason,
                            void *monitordestroy);

//###################################################################
/**Initializes Piecewise Linear FEM for diffusion solver.*/
int chi_diffusion::Solver::InitializePWLD(bool verbose)
{
  //Right now I am only doing one region at a time.
  //Later I want to support multiple regions with interfaces.
  chi_mesh::Region*     aregion = this->regions.back();
  grid = aregion->volume_mesh_continua.back();

  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  mesher = mesh_handler->volume_mesher;

  //================================================== Add pwl fem views
  if (verbose)
    chi_log.Log(LOG_0) << "Computing cell matrices";
  pwl_sdm = ((SpatialDiscretization_PWL*)(this->discretization));
  pwl_sdm->AddViewOfLocalContinuum(grid);
  pwl_sdm->AddViewOfNeighborContinuums(grid);
  MPI_Barrier(MPI_COMM_WORLD);

  //================================================== Reorder nodes
  if (verbose)
    chi_log.Log(LOG_0) << "Computing nodal reorderings for PWLD";
  ChiTimer t_reorder; t_reorder.Reset();
  ReorderNodesPWLD();
  auto domain_ownership = pwl_sdm->OrderNodesDFEM(grid);
  local_dof_count = domain_ownership.first;
  global_dof_count   = domain_ownership.second;

  MPI_Barrier(MPI_COMM_WORLD);
  if (verbose)
    chi_log.Log(LOG_0) << "Time taken during nodal reordering "
                       << t_reorder.GetTime()/1000.0;


  //================================================== Initialize field function
  //                                                   if empty
  pwld_phi_local.resize(local_dof_count);
  if (field_functions.empty())
  {
    auto initial_field_function = new chi_physics::FieldFunction(
      std::string("phi"),                           //Text name
      chi_physics_handler.fieldfunc_stack.size(),   //FF-id
      chi_physics::FieldFunctionType::DFEM_PWL,     //Type
      grid,                                         //Grid
      discretization,                               //Spatial Discretization
      1,                                            //Number of components
      1,                                            //Number of sets
      0,0,                                          //Ref component, ref set
      &pwld_cell_dof_array_address,                 //Dof block address
      &pwld_phi_local);                             //Data vector

    field_functions.push_back(initial_field_function);
    chi_physics_handler.fieldfunc_stack.push_back(initial_field_function);
  }

  //================================================== Setup timers
  if (verbose)
    chi_log.Log(LOG_0) << "Determining nodal connections";
  ChiTimer t_connect; t_connect.Reset();


  //================================================== Initialize nodal DOF
  //                                                   and connection info
  nodal_nnz_in_diag.resize(local_dof_count, 0);
  nodal_nnz_off_diag.resize(local_dof_count, 0);
  nodal_boundary_numbers.resize(grid->vertices.size(), 0);


  //================================================== Building sparsity pattern
  chi_log.Log(LOG_0) << "Building sparsity pattern.";
//  PWLDBuildSparsityPattern();
  pwl_sdm->BuildDFEMSparsityPattern(grid,
                                    nodal_nnz_in_diag,
                                    nodal_nnz_off_diag,
                                    domain_ownership);


  //================================================== Initialize x and b
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetSizes(x, local_dof_count, global_dof_count);CHKERRQ(ierr);
  ierr = VecSetType(x,VECMPI);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

  VecSet(x,0.0);
  VecSet(b,0.0);

  //################################################## Create matrix
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A, local_dof_count , local_dof_count,
                        global_dof_count, global_dof_count);CHKERRQ(ierr);
  ierr = MatSetType(A,MATMPIAIJ);CHKERRQ(ierr);

  //================================================== Allocate matrix memory
  chi_log.Log(LOG_0) << "Setting matrix preallocation.";
  MatMPIAIJSetPreallocation(A,0,nodal_nnz_in_diag.data(),
                              0,nodal_nnz_off_diag.data());
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
  MatSetUp(A);

  //================================================== Set up solver
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
  ierr = KSPSetOperators(ksp,A,A);
  ierr = KSPSetType(ksp,KSPCG);

  ierr = KSPGetPC(ksp,&pc);
  PCSetType(pc,PCHYPRE);

  PCHYPRESetType(pc,"boomeramg");

  //================================================== Setting Hypre parameters
  //The default HYPRE parameters used for polyhedra
  //seemed to have caused a lot of trouble for Slab
  //geometries. This section makes some custom options
  //per cell type
  auto first_cell = &grid->local_cells[0];

  if (first_cell->Type() == chi_mesh::CellType::SLAB)
  {
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_agg_nl 1");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_P_max 4");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_grid_sweeps_coarse 1");

    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_grid_sweeps_coarse 1");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_max_levels 25");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_coarsen_type HMIS");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_interp_type ext+i");

    PetscOptionsInsertString(NULL,"-options_left");
  }
  if (first_cell->Type() == chi_mesh::CellType::POLYGON)
  {
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_strong_threshold 0.6");

    //PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_agg_nl 1");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_P_max 4");

    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_grid_sweeps_coarse 1");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_max_levels 25");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_coarsen_type HMIS");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_interp_type ext+i");

    PetscOptionsInsertString(NULL,"-options_left");
  }
  if (first_cell->Type() == chi_mesh::CellType::POLYHEDRON)
  {
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_strong_threshold 0.8");

    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_agg_nl 1");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_P_max 4");

    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_grid_sweeps_coarse 1");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_max_levels 25");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_coarsen_type HMIS");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_interp_type ext+i");

    PetscOptionsInsertString(NULL,"-options_left");
  }
  PetscOptionsInsertString(NULL,options_string.c_str());
  PCSetFromOptions(pc);

  //=================================== Set up monitor
  if (verbose)
    ierr = KSPMonitorSet(ksp,&chi_diffusion::KSPMonitorAChiTech,NULL,NULL);

  KSPSetConvergenceTest(ksp,&DiffusionConvergenceTestNPT,NULL,NULL);

  //=================================== Setup verbose_info viewer
  if (chi_log.GetVerbosity()>= LOG_0VERBOSE_2)
    KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

  ierr = KSPSetTolerances(ksp,1.e-50,residual_tolerance,1.0e50,max_iters);
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

  return 0;
}