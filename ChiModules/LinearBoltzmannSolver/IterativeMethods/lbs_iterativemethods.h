#ifndef LBS_ITERATIVE_METHODS_H
#define LBS_ITERATIVE_METHODS_H

namespace lbs
{
  enum class IterativeMethod : int
  {
    NONE                     = 0,
    CLASSICRICHARDSON        = 1, ///< Otherwise known as Source Iteration
    CLASSICRICHARDSON_CYCLES = 2, ///< Source Iteration with Cycles support
    GMRES                    = 3, ///< GMRES iterative algorithm
    GMRES_CYCLES             = 4, ///< GMRES with Cycles support
    KRYLOV_RICHARDSON        = 5, ///< Richardson iteration
    KRYLOV_RICHARDSON_CYCLES = 6, ///< Richardson iteration with cycles support
    KRYLOV_GMRES             = 7, ///< GMRES iterative algorithm
    KRYLOV_GMRES_CYCLES      = 8, ///< GMRES with Cycles support
    KRYLOV_BICGSTAB          = 9, ///< BiCGStab iterative algorithm
    KRYLOV_BICGSTAB_CYCLES   = 10,///< BiCGStab with Cycles support
  };
}

#endif