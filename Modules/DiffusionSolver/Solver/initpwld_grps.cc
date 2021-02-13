#include "diffusion_solver.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiTimer/chi_timer.h>

#include <chi_mpi.h>
#include <chi_log.h>
#include <ChiPhysics/chi_physics.h>

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;
extern ChiPhysics&  chi_physics_handler;

PetscErrorCode
DiffusionConvergenceTestNPT(KSP ksp, PetscInt n, PetscReal rnorm,
                            KSPConvergedReason* convergedReason,
                            void *monitordestroy);

//###################################################################
/**Initializes Piecewise Linear FEM for diffusion solver.*/
int chi_diffusion::Solver::InitializePWLDGroups(bool verbose)
{
  //================================================== Add pwl fem views
  if (verbose)
    chi_log.Log(LOG_0) << "Computing cell matrices";
  pwl_sdm = std::static_pointer_cast<SpatialDiscretization_PWL>(this->discretization);
  pwl_sdm->PreComputeCellSDValues(grid);
  MPI_Barrier(MPI_COMM_WORLD);

  //================================================== Reorder nodes
  if (verbose)
    chi_log.Log(LOG_0) << "Computing nodal reorderings for PWLD";
  ChiTimer t_reorder; t_reorder.Reset();
  ReorderNodesPWLD();

  MPI_Barrier(MPI_COMM_WORLD);
  if (verbose)
    chi_log.Log(LOG_0) << "Time taken during nodal reordering "
                       << t_reorder.GetTime()/1000.0;

  //================================================== Initialize unknown manager
  unknown_manager.AddUnknown(chi_math::UnknownType::VECTOR_N, G);

  //================================================== Initialize field function
  //                                                   if empty
  pwld_phi_local.resize(local_dof_count * G);
  if (field_functions.empty())
  {
//    auto initial_field_function = new chi_physics::FieldFunction(
//      std::string("phi0"),                          //Text name
//      chi_physics_handler.fieldfunc_stack.size(),   //FF-id
//      chi_physics::FieldFunctionType::DFEM_PWL,     //Type
//      grid,                                         //Grid
//      discretization,                               //Spatial Discretization
//      1,                                            //Number of components
//      1,                                            //Number of sets
//      0,0,                                          //Ref component, ref set
//      &pwld_cell_dof_array_address,                 //Dof block address
//      &pwld_phi_local);                             //Data vector

    auto initial_field_function = new chi_physics::FieldFunction(
      std::string("phi"),                           //Text name
      discretization,                               //Spatial Discretization
      &pwld_phi_local,                              //Data vector
      unknown_manager);                             //Unknown Manager

    field_functions.push_back(initial_field_function);
    chi_physics_handler.fieldfunc_stack.push_back(initial_field_function);
  }

  //================================================== Setup timer
  if (verbose)
    chi_log.Log(LOG_0) << "Determining nodal connections";
  ChiTimer t_connect; t_connect.Reset();
  double t0 = 0.0;


  //================================================== Initialize nodal DOF
  //                                                   and connection info
  nodal_nnz_in_diag.resize(local_dof_count, 0);
  nodal_nnz_off_diag.resize(local_dof_count, 0);
  nodal_boundary_numbers.resize(grid->vertices.size(), 0);
  int total_nnz = 0;







  //================================================== First pass store pwld view
  chi_log.Log(LOG_0) << "Building sparsity pattern.";
  PWLDBuildSparsityPattern();


  //################################################## Initialize groupwise
  //                                                   Petsc data
  xg = new Vec[G];
  bg = new Vec[G];
  Ag = new Mat[G];
  kspg = new KSP[G];
  pcg = new PC[G];

  for (int gr=0; gr<G; gr++)
  {
    //=========================================== Initialize xg[gr] and bg[gr]
    ierr = VecCreate(PETSC_COMM_WORLD,&xg[gr]);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) xg[gr], "Solution");CHKERRQ(ierr);
    ierr = VecSetSizes(xg[gr], local_dof_count,
                       global_dof_count);CHKERRQ(ierr);
    ierr = VecSetType(xg[gr],VECMPI);CHKERRQ(ierr);
    ierr = VecDuplicate(xg[gr],&bg[gr]);CHKERRQ(ierr);

    //VecSet(xg[gr],0.0);
    VecSet(bg[gr],0.0);

    //=========================================== Create matrix
    ierr = MatCreate(PETSC_COMM_WORLD,&Ag[gr]);CHKERRQ(ierr);
    ierr = MatSetSizes(Ag[gr], local_dof_count,
                       local_dof_count,
                       global_dof_count, global_dof_count);CHKERRQ(ierr);
    ierr = MatSetType(Ag[gr],MATMPIAIJ);CHKERRQ(ierr);


    //=========================================== Allocate matrix memory
    MatMPIAIJSetPreallocation(Ag[gr],0,nodal_nnz_in_diag.data(),
                              0,nodal_nnz_off_diag.data());
    MatSetOption(Ag[gr], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetOption(Ag[gr],MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);
    MatSetUp(Ag[gr]);

    //================================================== Set up solver
    KSPCreate(PETSC_COMM_WORLD,&kspg[gr]);
    KSPSetOperators(kspg[gr],Ag[gr],Ag[gr]);
    KSPSetType(kspg[gr],KSPCG);

    KSPGetPC(kspg[gr],&pcg[gr]);
    PCSetType(pcg[gr],PCHYPRE);

    PCHYPRESetType(pcg[gr],"boomeramg");

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
    PCSetFromOptions(pcg[gr]);


    //=================================== Set up monitor
    if (verbose)
      ierr = KSPMonitorSet(kspg[gr],&chi_diffusion::KSPMonitorAChiTech,NULL,NULL);

    KSPSetConvergenceTest(kspg[gr],&DiffusionConvergenceTestNPT,NULL,NULL);

    //=================================== Setup verbose_info viewer
    if (chi_log.GetVerbosity()>= LOG_0VERBOSE_2)
      KSPView(kspg[gr],PETSC_VIEWER_STDOUT_WORLD);

    ierr = KSPSetTolerances(kspg[gr],1.e-50,residual_tolerance,1.0e50,max_iters);
    ierr = KSPSetInitialGuessNonzero(kspg[gr],PETSC_TRUE);
  }


  return 0;
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

  if (relative_residual < tol)
    *convergedReason = KSP_CONVERGED_RTOL;

  return KSP_CONVERGED_ITERATING;
}