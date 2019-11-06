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
int chi_diffusion::Solver::InitializePWLDGroups(bool verbose)
{
  //Right now I am only doing one region at a time.
  //Later I want to support multiple regions with interfaces.
  chi_mesh::Region*     aregion = this->regions.back();
  grid = aregion->volume_mesh_continua.back();

  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  mesher = mesh_handler->volume_mesher;

  int num_nodes = grid->nodes.size();

  //================================================== Reorder nodes
  if (verbose)
    chi_log.Log(LOG_0) << "Computing nodal reorderings for PWLD";
  ChiTimer t_reorder; t_reorder.Reset();
  ReorderNodesPWLD();

  MPI_Barrier(MPI_COMM_WORLD);
  if (verbose)
    chi_log.Log(LOG_0) << "Time taken during nodal reordering "
                       << t_reorder.GetTime()/1000.0;


  //================================================== Initialize field function
  //                                                   if empty
  pwld_phi_local.resize(pwld_local_dof_count*G);
  if (field_functions.size() == 0)
  {
    chi_physics::FieldFunction* initial_field_function =
      new chi_physics::FieldFunction;
    initial_field_function->text_name = std::string("phi0");
    initial_field_function->grid = grid;
    initial_field_function->spatial_discretization = discretization;
    initial_field_function->id = chi_physics_handler.fieldfunc_stack.size();

    initial_field_function->type = FF_SDM_PWLD;
    initial_field_function->num_grps = 1;
    initial_field_function->num_moms = 1;
    initial_field_function->grp = 0;
    initial_field_function->mom = 0;
    initial_field_function->field_vector_local = &pwld_phi_local;
    initial_field_function->local_cell_dof_array_address =
      &pwld_cell_dof_array_address;

    field_functions.push_back(initial_field_function);
    chi_physics_handler.fieldfunc_stack.push_back(initial_field_function);
  }
  else
  {
    size_t num_ff = field_functions.size();
    for (int ff=0; ff<num_ff; ff++)
    {
      chi_physics::FieldFunction* cur_ff = field_functions[ff];
      cur_ff->grid                   = grid;
      cur_ff->spatial_discretization = discretization;

      cur_ff->type = FF_SDM_PWLD;
      cur_ff->num_grps = 1;
      cur_ff->num_moms = 1;
      cur_ff->grp = 0;
      cur_ff->mom = 0;
      cur_ff->field_vector_local = &pwld_phi_local;
      cur_ff->local_cell_dof_array_address = &pwld_cell_dof_array_address;
    }
  }

  //================================================== Setup timer
  if (verbose)
    chi_log.Log(LOG_0) << "Determining nodal connections";
  ChiTimer t_connect; t_connect.Reset();
  double t0 = 0.0;


  //================================================== Initialize nodal DOF
  //                                                   and connection info
  nodal_nnz_in_diag.resize(pwld_local_dof_count,0);
  nodal_nnz_off_diag.resize(pwld_local_dof_count,0);
  nodal_boundary_numbers.resize(grid->nodes.size(),0);
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
    ierr = VecSetSizes(xg[gr],pwld_local_dof_count,
                       pwld_global_dof_count);CHKERRQ(ierr);
    ierr = VecSetType(xg[gr],VECMPI);CHKERRQ(ierr);
    ierr = VecDuplicate(xg[gr],&bg[gr]);CHKERRQ(ierr);

    //VecSet(xg[gr],0.0);
    VecSet(bg[gr],0.0);

    //=========================================== Create matrix
    ierr = MatCreate(PETSC_COMM_WORLD,&Ag[gr]);CHKERRQ(ierr);
    ierr = MatSetSizes(Ag[gr],pwld_local_dof_count,
                       pwld_local_dof_count,
                       pwld_global_dof_count,pwld_global_dof_count);CHKERRQ(ierr);
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
    int first_cell_g_index = grid->local_cell_glob_indices[0];
    auto first_cell = grid->cells[first_cell_g_index];

    if (first_cell->Type() == chi_mesh::CellType::CELL_NEWBASE)
    {
      auto cell_base = dynamic_cast<chi_mesh::CellBase*>(first_cell);

      if (cell_base->Type2() == chi_mesh::CellType::SLABV2)
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
      if (cell_base->Type2() == chi_mesh::CellType::POLYGONV2)
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
      if (cell_base->Type2() == chi_mesh::CellType::POLYHEDRONV2)
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


    }
    PetscOptionsInsertString(NULL,options_string.c_str());
    PCSetFromOptions(pcg[gr]);


    //=================================== Set up monitor
    if (verbose)
      ierr = KSPMonitorSet(kspg[gr],&chi_diffusion::KSPMonitorAChiTech,NULL,NULL);

    KSPSetConvergenceTest(kspg[gr],&DiffusionConvergenceTestNPT,NULL,NULL);

    //=================================== Setup verbose viewer
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