#include "diffusion_solver.h"

#include "ChiPhysics/FieldFunction2/fieldfunction2.h"

#include "chi_runtime.h"

#include "chi_log.h"
#include "chi_mpi.h"

#include "ChiTimer/chi_timer.h"

//###################################################################
/**Initializes the diffusion solver using the PETSc library.*/
int chi_diffusion::Solver::Initialize(bool verbose)
{
  chi::log.Log() << "\n"
                     << chi::program_timer.GetTimeString() << " "
                     << TextName() << ": Initializing Diffusion solver ";
  this->verbose_info = verbose;

  if (not common_items_initialized)
    InitializeCommonItems(); //Mostly boundaries

  chi_objects::ChiTimer t_init; t_init.Reset();

  auto sdm_string = basic_options("discretization_method").StringValue();
  {
    using namespace chi_math::finite_element;
    if      (sdm_string == "PWLC")
    {
      discretization =
        chi_math::SpatialDiscretization_PWLC::New(grid, COMPUTE_UNIT_INTEGRALS);
      unknown_manager.AddUnknown(chi_math::UnknownType::SCALAR);
    }
    else if (sdm_string == "PWLD_MIP")
    {
      discretization =
        chi_math::SpatialDiscretization_PWLD::New(grid, COMPUTE_UNIT_INTEGRALS);
      unknown_manager.AddUnknown(chi_math::UnknownType::SCALAR);
    }
//    else if (sdm_string == "PWLD_MIP_GAGG")
//    {
//      discretization =
//        chi_math::SpatialDiscretization_PWLD::New(grid, COMPUTE_UNIT_INTEGRALS);
//      unknown_manager.AddUnknown(chi_math::UnknownType::VECTOR_N, G);
//    }
    else
      throw std::invalid_argument(
        TextName() + ": Invalid spatial discretization method, " +
        sdm_string + ", specified.");
  }

  MPI_Barrier(MPI_COMM_WORLD);
  auto& sdm = discretization;

  //============================================= Get DOF counts
  local_dof_count = sdm->GetNumLocalDOFs(unknown_manager);
  global_dof_count = sdm->GetNumGlobalDOFs(unknown_manager);
  chi::log.Log()
    << TextName() << ": Global number of DOFs="
    << global_dof_count;


  //================================================== Initialize discretization
  //                                                   method
  if (field_functions2.empty())
  {
    auto& sdm_ptr = discretization;
    std::string solver_name;
    if (not TextName().empty()) solver_name = TextName() + "-";

    std::string text_name = solver_name + "phi";

    using namespace chi_math;
    auto initial_field_function =
      std::make_shared<chi_physics::FieldFunction2>(
        text_name,                     //Text name
        sdm_ptr,                       //Spatial Discretization
        Unknown(UnknownType::SCALAR)); //Unknown/Variable

    field_functions2.push_back(initial_field_function);
    chi::fieldfunc2_stack.push_back(initial_field_function);
  }//if not ff set


  //================================================== Determine nodal DOF
  chi::log.Log() << "Building sparsity pattern.";
  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm->BuildSparsityPattern(nodal_nnz_in_diag,
                            nodal_nnz_off_diag,
                            unknown_manager);

  chi::log.Log()
    << chi::program_timer.GetTimeString() << " "
    << TextName() << ": Diffusion Solver initialization time "
    << t_init.GetTime()/1000.0 << std::endl;

  //================================================== Initialize x and b
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetSizes(x, static_cast<PetscInt>(local_dof_count),
                        static_cast<PetscInt>(global_dof_count));CHKERRQ(ierr);
  ierr = VecSetType(x,VECMPI);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

  VecSet(x,0.0);
  VecSet(b,0.0);

  //################################################## Create matrix
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A, static_cast<PetscInt>(local_dof_count),
                        static_cast<PetscInt>(local_dof_count),
                        static_cast<PetscInt>(global_dof_count),
                        static_cast<PetscInt>(global_dof_count));CHKERRQ(ierr);
  ierr = MatSetType(A,MATMPIAIJ);CHKERRQ(ierr);

  //================================================== Allocate matrix memory
  chi::log.Log() << "Setting matrix preallocation.";
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

  KSPSetConvergenceTest(ksp,&chi_diffusion::DiffusionConvergenceTestNPT,NULL,NULL);

  ierr = KSPSetTolerances(ksp,
                          1.e-50,
                          basic_options("residual_tolerance").FloatValue(),
                          1.0e50,
                          basic_options("max_iters").IntegerValue());
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

  return false;
}