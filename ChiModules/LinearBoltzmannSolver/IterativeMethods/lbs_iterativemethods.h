#ifndef LBS_ITERATIVE_METHODS_H
#define LBS_ITERATIVE_METHODS_H

namespace LinearBoltzmann
{
  enum class IterativeMethod : int
  {
    NONE                     = 0,
    CLASSICRICHARDSON        = 1,
    CLASSICRICHARDSON_CYCLES = 2,
    GMRES                    = 3,
    GMRES_CYCLES             = 4
  };
}

#endif