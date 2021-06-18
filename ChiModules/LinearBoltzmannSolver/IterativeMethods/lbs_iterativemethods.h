#ifndef LBS_ITERATIVE_METHODS_H
#define LBS_ITERATIVE_METHODS_H

namespace LinearBoltzmann
{
  enum class IterativeMethod : int
  {
    NONE                     = 0,
    CLASSICRICHARDSON        = 1, ///< Otherwise known as Source Iteration
    CLASSICRICHARDSON_CYCLES = 2, ///< Source Iteration with Cycles support
    GMRES                    = 3, ///< GMRES iterative algorithm
    GMRES_CYCLES             = 4  ///< GMRES with Cycles support
  };
}

#endif