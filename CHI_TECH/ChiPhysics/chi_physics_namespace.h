#ifndef _chi_physics_namespace_h
#define _chi_physics_namespace_h
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
    EXISTING     = 23
  };

  class FieldFunction;
  class Solver;
  class Material;
  class MaterialProperty;

  class ScalarValue;
  class TransportCrossSections;
  class IsotropicMultiGrpSource;

  //03 Utils
  std::string GetPETScConvergedReasonstring(KSPConvergedReason reason);
}


#endif