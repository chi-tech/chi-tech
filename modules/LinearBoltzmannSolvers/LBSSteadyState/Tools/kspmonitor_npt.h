#ifndef LBS_KSP_MONITORS_H
#define LBS_KSP_MONITORS_H

#include <petscksp.h>

namespace lbs
{
PetscErrorCode KSPConvergenceTest(
  KSP ksp, PetscInt n, PetscReal rnorm,
  KSPConvergedReason* convergedReason, void* monitordestroy);
}

#endif //LBS_KSP_MONITORS_H