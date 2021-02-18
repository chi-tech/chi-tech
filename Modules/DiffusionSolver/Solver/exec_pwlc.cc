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
/**Builds the matrix using the PWLC discretization method.*/
int chi_diffusion::Solver::ExecutePWLC(bool suppress_assembly,
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
                       << solver_name << ": Assembling A locally";

  //================================================== Loop over locally owned
  //                                                   cells
  if (!suppress_assembly)
    for (auto& cell : grid->local_cells)
      CFEM_Assemble_A_and_b(&cell, gi);




  if (!suppress_assembly)
    chi_log.Log(LOG_0) << chi_program_timer.GetTimeString() << " "
                       << solver_name << ": Done Assembling A locally";
  MPI_Barrier(MPI_COMM_WORLD);

  //=================================== Call matrix assembly
  if (verbose_info || chi_log.GetVerbosity() >= LOG_0VERBOSE_1)
    chi_log.Log(LOG_0)
      << chi_program_timer.GetTimeString() << " "
      << solver_name << ": Communicating matrix assembly";

  if (!suppress_assembly)
  {
    chi_log.Log(LOG_0) << chi_program_timer.GetTimeString() << " "
                       << solver_name << ": Assembling A globally";
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
                       << solver_name << ": Diagonal check";
    PetscBool missing_diagonal;
    PetscInt  row;
    MatMissingDiagonal(A,&missing_diagonal,&row);
    if (missing_diagonal)
      chi_log.Log(LOG_ALLERROR) << chi_program_timer.GetTimeString() << " "
                                << solver_name << ": Missing diagonal detected";

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
                       << solver_name << ": Assembling x and b";
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
      << solver_name
      << ": Setting up solver and preconditioner\n";
    PCSetUp(pc);
    KSPSetUp(ksp);
  }
  else
  {
    if (verbose_info || chi_log.GetVerbosity() >= LOG_0VERBOSE_1)
      chi_log.Log(LOG_0)
        << chi_program_timer.GetTimeString() << " "
        << solver_name << ": Solving system\n";
    t_solve.Reset();
    PCSetUp(pc);
    KSPSetUp(ksp);
    KSPSolve(ksp,b,x);
    time_solve = t_solve.GetTime()/1000.0;










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
      int its;
      ierr = KSPGetIterationNumber(ksp,&its);
      chi_log.Log(LOG_0)
        << chi_program_timer.GetTimeString() << " "
        << solver_name
        << ": Number of iterations =" << its;

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
