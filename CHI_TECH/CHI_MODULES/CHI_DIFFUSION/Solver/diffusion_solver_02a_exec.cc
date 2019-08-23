#include "diffusion_solver.h"

#include <chi_log.h>
extern CHI_LOG chi_log;

//###################################################################
/**Executes the diffusion solver using the PETSc library.*/
int chi_diffusion::Solver::ExecuteS(bool suppress_assembly,
                                    bool suppress_solve)
{
  if (fem_method == PWLC)
    ExecutePWLC(suppress_assembly,
                suppress_solve);
  else if (fem_method == PWLD_MIP)
    ExecutePWLD_MIP(suppress_assembly,
                    suppress_solve);
  else if (fem_method == PWLD_MIP_GRPS)
    ExecutePWLD_MIP_GRPS(suppress_assembly,
                         suppress_solve);
  else if (fem_method == PWLD_MIP_GAGG)
    ExecutePWLD_MIP_GAGG(suppress_assembly,
                         suppress_solve);
  else
  {
    chi_log.Log(LOG_0)
      << "Diffusion Solver: Finite Element Discretization "
         "method not specified.";
    exit(EXIT_FAILURE);
  }

  return 0;
}