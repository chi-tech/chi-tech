#include "diffusion_solver.h"

#include "ChiTimer/chi_timer.h"

#include "chi_mpi.h"
#include "chi_log.h"
#include "ChiPhysics/chi_physics.h"

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;
extern ChiPhysics&  chi_physics_handler;

extern ChiTimer chi_program_timer;

//###################################################################
/**Executes the diffusion solver using the PETSc library.*/
int chi_diffusion::Solver::ExecuteS(bool suppress_assembly,
                                    bool suppress_solve)
{
  t_assembly.Reset();

  if (chi_physics_handler.material_stack.empty())
  {
    chi_log.Log(LOG_0ERROR)
      << "No materials added to simulation. Add materials.";
    exit(0);
  }

  VecSet(x,0.0);
  VecSet(b,0.0);

  if (!suppress_assembly)
    chi_log.Log(LOG_0) << chi_program_timer.GetTimeString() << " "
                       << TextName() << ": Assembling A locally";

  //================================================== Loop over locally owned
  //                                                   cells
  if (fem_method == PWLC)
  {
    if (!suppress_assembly)
      for (auto& cell : grid->local_cells)
        CFEM_Assemble_A_and_b(cell, gi);
    else {}
  }
  else if (fem_method == PWLD_MIP)
  {
    if (!suppress_assembly)
      for (auto& cell : grid->local_cells)
        PWLD_Assemble_A_and_b(cell,gi);
    else
      for (auto& cell : grid->local_cells)
        PWLD_Assemble_b(cell,gi);
  }
  else if (fem_method == PWLD_MIP_GAGG)
  {
    if (!suppress_assembly)
      for (auto& cell : grid->local_cells)
        PWLD_Assemble_A_and_b_GAGG(cell);
    else
      for (auto& cell : grid->local_cells)
        PWLD_Assemble_b_GAGG(cell);
  }
  else
  {
    chi_log.Log(LOG_0)
      << "Diffusion Solver: Finite Element Discretization "
         "method not specified.";
    exit(EXIT_FAILURE);
  }

  if (!suppress_assembly)
    chi_log.Log(LOG_0) << chi_program_timer.GetTimeString() << " "
                       << TextName() << ": Done Assembling A locally";
  MPI_Barrier(MPI_COMM_WORLD);

  //=================================== Call matrix assembly
  if (verbose_info || chi_log.GetVerbosity() >= LOG_0VERBOSE_1)
    chi_log.Log(LOG_0)
      << chi_program_timer.GetTimeString() << " "
      << TextName() << ": Communicating matrix assembly";

  if (!suppress_assembly)
  {
    chi_log.Log(LOG_0) << chi_program_timer.GetTimeString() << " "
                       << TextName() << ": Assembling A globally";
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    //================================= Matrix symmetry check
//    PetscBool is_symmetric;
//    ierr = MatIsSymmetric(A,1.0e-4,&is_symmetric);
//    if (!is_symmetric)
//    {
//      chi_log.Log(LOG_0WARNING)
//        << "Assembled matrix is not symmetric";
//    }

    //================================= Matrix diagonal check
    chi_log.Log(LOG_0) << chi_program_timer.GetTimeString() << " "
                       << TextName() << ": Diagonal check";
    PetscBool missing_diagonal;
    PetscInt  row;
    MatMissingDiagonal(A,&missing_diagonal,&row);
    if (missing_diagonal)
      chi_log.Log(LOG_ALLERROR) << chi_program_timer.GetTimeString() << " "
                                << TextName() << ": Missing diagonal detected";

    //================================= Matrix sparsity info
    MatInfo info;
    ierr = MatGetInfo(A,MAT_GLOBAL_SUM,&info);

    chi_log.Log(LOG_0) << "Number of mallocs used = " << info.mallocs
                       << "\nNumber of non-zeros allocated = "
                       << info.nz_allocated
                       << "\nNumber of non-zeros used = "
                       << info.nz_used
                       << "\nNumber of unneeded non-zeros = "
                       << info.nz_unneeded;
  }
  if (verbose_info || chi_log.GetVerbosity() >= LOG_0VERBOSE_1)
    chi_log.Log(LOG_0) << chi_program_timer.GetTimeString() << " "
                       << TextName() << ": Assembling x and b";
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  time_assembly = t_assembly.GetTime()/1000.0;

  //=================================== Execute solve
  if (suppress_solve)
  {
    chi_log.Log(LOG_0)
      << chi_program_timer.GetTimeString() << " "
      << TextName()
      << ": Setting up solver and preconditioner\n";
    PCSetUp(pc);
    KSPSetUp(ksp);
  }
  else
  {
    if (verbose_info || chi_log.GetVerbosity() >= LOG_0VERBOSE_1)
      chi_log.Log(LOG_0)
        << chi_program_timer.GetTimeString() << " "
        << TextName() << ": Solving system\n";
    t_solve.Reset();
    PCSetUp(pc);
    KSPSetUp(ksp);
    KSPSolve(ksp,b,x);
    time_solve = t_solve.GetTime()/1000.0;

    //=================================== Populate field vector
    if (fem_method == PWLD_MIP or fem_method == PWLD_MIP_GAGG)
    {
      const double* x_ref;
      VecGetArrayRead(x,&x_ref);

      for (int i=0; i < local_dof_count; i++)
        pwld_phi_local[i] = x_ref[i];

      VecRestoreArrayRead(x,&x_ref);
    }

    //=================================== Get convergence reason
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp,&reason);
    if (verbose_info || reason != KSP_CONVERGED_RTOL)
      chi_log.Log(LOG_0)
        << "Convergence reason: "
        << chi_physics::GetPETScConvergedReasonstring(reason);

    //=================================== Location wise view
    if (chi_mpi.location_id == 0)
    {
      int64_t its;
      ierr = KSPGetIterationNumber(ksp,&its);
      chi_log.Log(LOG_0)
        << chi_program_timer.GetTimeString() << " "
        << TextName() << "[g=" << gi << "-" << gi+G-1
        << "]: Number of iterations =" << its;

      if (verbose_info || chi_log.GetVerbosity() >= LOG_0VERBOSE_1)
      {
        chi_log.Log(LOG_0) << "Timing:";
        chi_log.Log(LOG_0) << "Assembling the matrix: " << time_assembly;
        chi_log.Log(LOG_0) << "Solving the system   : " << time_solve;
      }
    }

    if (verbose_info || chi_log.GetVerbosity() >= LOG_0VERBOSE_1)
      chi_log.Log(LOG_0) << "Diffusion Solver execution completed!\n";
  }//if not suppressed solve

  return 0;
}