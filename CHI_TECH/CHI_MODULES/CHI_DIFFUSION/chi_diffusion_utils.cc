#include "chi_diffusion.h"

#include "../../CHI_MPI/chi_mpi.h"
#include <chi_log.h>

extern CHI_MPI chi_mpi;
extern CHI_LOG chi_log;

//###################################################################
/**Customized monitor for PETSc Krylov sub-space solvers.*/
PetscErrorCode chi_diffusion::
KSPMonitorAChiTech(KSP ksp, PetscInt n, PetscReal rnorm, void *monitordestroy)
{

  Vec Rhs;
  KSPGetRhs(ksp,&Rhs);
  double rhs_norm;
  VecNorm(Rhs,NORM_2,&rhs_norm);
  if (rhs_norm < 1.0e-25)
    rhs_norm = 1.0;

  if (chi_mpi.location_id == 0)
  {

    char buff[100];
    if (rnorm/rhs_norm < 1.0e-2)
    {
      snprintf(buff,100,"Diffusion iteration %4d - Residual %.3e\n",n,rnorm/rhs_norm);
    }
    else
    {
      snprintf(buff,100,"Diffusion iteration %4d - Residual %.7f\n",n,rnorm/rhs_norm);
    }

    chi_log.Log(LOG_0) << buff;
  }



  return 0;
}