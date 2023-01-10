#ifndef CHI_PHYSICS_NAMESPACE_H
#define CHI_PHYSICS_NAMESPACE_H
/**\defgroup LuaPhysics C Physics*/

#include <petscksp.h>

namespace chi_physics
{
  enum class OperationType
  {
    SINGLE_VALUE = 0,
    FROM_ARRAY   = 1,
    SIMPLEXS0    = 20,
    SIMPLEXS1    = 21,
    PDT_XSFILE   = 22,
    EXISTING     = 23,
    CHI_XSFILE   = 24
  };

  class FieldFunction;
  class Solver;

  //03 Utils
  std::string GetPETScConvergedReasonstring(KSPConvergedReason reason);
}


#endif