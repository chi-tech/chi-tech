#include <petscksp.h>

PetscErrorCode KSPMonitorNPT(KSP ksp,
                             PetscInt n,
                             PetscReal rnorm,
                             void *monitordestroy);

PetscErrorCode KSPConvergenceTestNPT(
                KSP ksp, PetscInt n, PetscReal rnorm,
                KSPConvergedReason* convergedReason, void *monitordestroy);