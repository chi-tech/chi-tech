#include "chi_diffusion.h"

#include <iomanip>

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiMPI/chi_mpi.h"
extern ChiMPI& chi_mpi;

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
    const auto ksp_name = "Diffusion";

    std::stringstream buff;
    buff
      << ksp_name
      << " iteration "
      << std::setw(4) << n
      << " - Residual "
      << std::scientific << std::setprecision(7) << rnorm / rhs_norm
      << std::endl;

    chi_log.Log(LOG_0) << buff.str();
  }



  return 0;
}