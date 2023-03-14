#ifndef CHITECH_LBS_K_RESIDUAL_FUNCTION_H
#define CHITECH_LBS_K_RESIDUAL_FUNCTION_H

#include <petscsnes.h>

namespace lbs
{
  PetscErrorCode SNESKResidualFunction(SNES, Vec, Vec, void *);
}

#endif //CHITECH_LBS_K_RESIDUAL_FUNCTION_H
