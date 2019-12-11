#ifndef _chi_physics_namespace_h
#define _chi_physics_namespace_h
/**\defgroup LuaPhysics C Physics*/



//Operation indices
#define SINGLE_VALUE              0
#define FROM_ARRAY                1
#define SIMPLEXS0                 20
#define SIMPLEXS1                 21
#define PDT_XSFILE                22
#define EXISTING                  23


namespace chi_physics
{
  class FieldFunction;
  class Solver;
  class Material;
  class MaterialProperty;

  class ThermalConductivity;
  class ScalarValue;
  class TransportCrossSections;
  class IsotropicMultiGrpSource;
}


#endif