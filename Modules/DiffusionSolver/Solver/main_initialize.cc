#include "diffusion_solver.h"



#include "ChiTimer/chi_timer.h"
#include "chi_mpi.h"
#include "chi_log.h"
#include "ChiPhysics/chi_physics.h"

extern ChiTimer chi_program_timer;
extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;
extern ChiPhysics&  chi_physics_handler;

PetscErrorCode
DiffusionConvergenceTestNPT(KSP ksp, PetscInt n, PetscReal rnorm,
                            KSPConvergedReason* convergedReason,
                            void *monitordestroy);

//###################################################################
/**Initializes the diffusion solver using the PETSc library.*/
int chi_diffusion::Solver::Initialize(bool verbose)
{
  chi_log.Log(LOG_0) << "\n"
                     << chi_program_timer.GetTimeString() << " "
                     << solver_name << ": Initializing Diffusion solver ";
  this->verbose_info = verbose;

  if (regions.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_diffusion::Solver::Initialize: No region added to solver.";
    exit(EXIT_FAILURE);
  }

  auto& region = regions.back();
  grid = region->GetGrid();



  if (not common_items_initialized)
    InitializeCommonItems(); //Mostly boundaries

  ChiTimer t_init; t_init.Reset();

  switch (fem_method)
  {
    using namespace chi_math::finite_element;
    case PWLC: {
      discretization =
        SpatialDiscretization_PWLC::New(grid, COMPUTE_UNIT_INTEGRALS);
      unknown_manager.AddUnknown(chi_math::UnknownType::SCALAR);
      break;
    }
    case PWLD_MIP: {
      discretization =
        SpatialDiscretization_PWL::New(grid, COMPUTE_UNIT_INTEGRALS);
      unknown_manager.AddUnknown(chi_math::UnknownType::SCALAR);
      break;
    }
    case PWLD_MIP_GAGG: {
      discretization =
        SpatialDiscretization_PWL::New(grid, COMPUTE_UNIT_INTEGRALS);
      unknown_manager.AddUnknown(chi_math::UnknownType::VECTOR_N, G);
      break;
    }
    default:
    {
      chi_log.Log(LOG_0)
        << "Diffusion Solver: Finite Element Discretization "
           "method not specified.";
      exit(EXIT_FAILURE);
    }
  }//switch fem_method
  MPI_Barrier(MPI_COMM_WORLD);
  auto& sdm = discretization;

  //============================================= Get DOF counts
  local_dof_count = sdm->GetNumLocalDOFs(grid, unknown_manager);
  global_dof_count = sdm->GetNumGlobalDOFs(grid, unknown_manager);
  chi_log.Log(LOG_0)
    << solver_name << ": Global number of DOFs="
    << global_dof_count;


  //================================================== Initialize discretization
  //                                                   method
  if (field_functions.empty())
  {
    switch (fem_method)
    {
      case PWLC:
      {
        auto initial_field_function = new chi_physics::FieldFunction(
          std::string("phi"),                           //Text name
          discretization,                               //Spatial Discretization
          &x,                                           //Data vector
          unknown_manager);                             //Unknown Manager

        field_functions.push_back(initial_field_function);
        chi_physics_handler.fieldfunc_stack.push_back(initial_field_function);
        break;
      }
      case PWLD_MIP:
      case PWLD_MIP_GAGG:
      {
        pwld_phi_local.resize(local_dof_count);
        if (field_functions.empty())
        {
          auto initial_field_function = new chi_physics::FieldFunction(
            std::string("phi"),                           //Text name
            discretization,                               //Spatial Discretization
            &pwld_phi_local,                              //Data vector
            unknown_manager);                             //Unknown Manager

          field_functions.push_back(initial_field_function);
          chi_physics_handler.fieldfunc_stack.push_back(initial_field_function);
        }
        break;
      }
    }//switch fem_method
  }//if not ff set


  //================================================== Determine nodal DOF
  chi_log.Log(LOG_0) << "Building sparsity pattern.";
  std::vector<int> nodal_nnz_in_diag;
  std::vector<int> nodal_nnz_off_diag;
  sdm->BuildSparsityPattern(grid,
                            nodal_nnz_in_diag,
                            nodal_nnz_off_diag,
                            unknown_manager);

  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString() << " "
    << solver_name << ": Diffusion Solver initialization time "
    << t_init.GetTime()/1000.0 << std::endl;

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
  ierr = MatSetSizes(A, local_dof_count, local_dof_count,
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

  //================================================== Set up preconditioner
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
  }
  PetscOptionsInsertString(NULL,options_string.c_str());
  PCSetFromOptions(pc);

  //=================================== Set up monitor
  if (verbose)
    ierr = KSPMonitorSet(ksp,&chi_diffusion::KSPMonitorAChiTech,NULL,NULL);

  KSPSetConvergenceTest(ksp,&DiffusionConvergenceTestNPT,NULL,NULL);

  ierr = KSPSetTolerances(ksp,1.e-50,residual_tolerance,1.0e50,max_iters);
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

  return false;
}

//###################################################################
/**Customized convergence test.*/
PetscErrorCode
DiffusionConvergenceTestNPT(KSP ksp, PetscInt n, PetscReal rnorm,
                            KSPConvergedReason* convergedReason, void *monitordestroy)
{
  //======================================================= Compute rhs norm
  Vec Rhs;
  KSPGetRhs(ksp,&Rhs);
  double rhs_norm;
  VecNorm(Rhs,NORM_2,&rhs_norm);
  if (rhs_norm < 1.0e-25)
    rhs_norm = 1.0;

  //======================================================= Compute test criterion
  double tol;
  int    maxIts;
  KSPGetTolerances(ksp,NULL,&tol,NULL,&maxIts);


  double relative_residual = rnorm/rhs_norm;

  chi_log.Log(LOG_0) << "Iteration " << n << " Residual " << rnorm/rhs_norm;

  if (relative_residual < tol)
    *convergedReason = KSP_CONVERGED_RTOL;

  return KSP_CONVERGED_ITERATING;
}