#include "diffusion_solver.h"

#include <PiecewiseLinear/CellViews/pwl_slab.h>
#include <PiecewiseLinear/CellViews/pwl_polygon.h>
#include <PiecewiseLinear/CellViews/pwl_polyhedron.h>
#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/chi_volumemesher.h>

#include "../Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "../Boundaries/chi_diffusion_bndry_reflecting.h"

#include <ChiTimer/chi_timer.h>

#include <chi_mpi.h>
#include <chi_log.h>
#include <ChiPhysics/chi_physics.h>
extern ChiMPI chi_mpi;
extern ChiLog chi_log;
extern ChiPhysics chi_physics_handler;


//###################################################################
/**Builds the matrix using the PWLC discretization method.*/
int chi_diffusion::Solver::ExecutePWLC(bool suppress_assembly,
                                       bool suppress_solve)
{
  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  mesher    = mesh_handler->volume_mesher;


  //################################################## Assemble Amatrix
  chi_log.Log(LOG_0) << "Diffusion Solver: Assembling A and b";
  chi_log.Log(LOG_0) << "Diffusion Solver: Local matrix instructions";
  t_assembly.Reset();

  std::vector<int> boundary_nodes,boundary_numbers;

  if (chi_physics_handler.material_stack.empty())
  {
    chi_log.Log(LOG_0ERROR)
      << "No materials added to simulation. Add materials.";
    exit(0);
  }

  //================================================== Setting references
  xref = x;
  bref = b;
  Aref = A;

  //================================================== Loop over locally owned
  //                                                   cells
  for (auto& cell : grid->local_cells)
  {
    if (!suppress_assembly)
      CFEM_Assemble_A_and_b(cell.cell_global_id, &cell, gi);
  }

  //=================================== Call matrix assembly
  chi_log.Log(LOG_0) << "Diffusion Solver: Communicating matrix assembly";

  if (!suppress_assembly)
  {
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    VecAssemblyBegin(x);
    VecAssemblyEnd(x);

    //================================= Matrix symmetry check
    PetscBool is_symmetric;
    ierr = MatIsSymmetric(A,1.0e-4,&is_symmetric);
    if (!is_symmetric)
    {
      chi_log.Log(LOG_0WARNING)
        << "Assembled matrix is not symmetric";
    }
  }
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  time_assembly = t_assembly.GetTime()/1000.0;

  //=================================== Execute solve
  if (suppress_solve)
  {
    chi_log.Log(LOG_0) << "Diffusion Solver: Setting up solver\n";
    PCSetUp(pc);
    KSPSetUp(ksp);
  }
  else
  {
    chi_log.Log(LOG_0) << "Diffusion Solver: Solving system\n";
    t_solve.Reset();
    PCSetUp(pc);
    KSPSetUp(ksp);
    KSPSolve(ksp,b,x);
    time_solve = t_solve.GetTime()/1000.0;

    //=================================== Get convergence reason
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp,&reason);
    chi_log.Log(LOG_0) << "Convergence reason: " << reason;


    //=================================== Location wise view
    if (chi_mpi.location_id == 0)
    {
      int its;
      ierr = KSPGetIterationNumber(ksp,&its);
      chi_log.Log(LOG_0) << "Diffusion Solver: Number of iterations =" << its;
      chi_log.Log(LOG_0) << "Timing:";
      chi_log.Log(LOG_0) << "Assembling the matrix: " << time_assembly;
      chi_log.Log(LOG_0) << "Solving the system   : " << time_solve;
    }

    chi_log.Log(LOG_0) << "Diffusion Solver execution completed!\n";
  }//if not suppressed solve

  return 0;
}