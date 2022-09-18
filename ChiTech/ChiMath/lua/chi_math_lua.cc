#include "chi_math_lua.h"

#include "ChiMath/Quadratures/LegendrePoly/lua/legendre_lua.h"
#include "ChiMath/Quadratures/lua/quadratures_lua.h"
#include "ChiMath/Quadratures/SLDFESQ/lua/sldfe_lua.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)
#define LUA_CMACRO1(x,y) \
        lua_pushnumber(L, y); \
        lua_setglobal(L, #x)

void chi_math::lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiLegendre);
  LUA_FMACRO1(chiLegendreDerivative);
  LUA_FMACRO1(chiYlm);

  LUA_FMACRO1(chiCreateCustomAngularQuadrature);
  LUA_FMACRO1(chiCreateCylindricalProductQuadrature);
  LUA_FMACRO1(chiCreateSphericalProductQuadrature);
  LUA_FMACRO1(chiCreateProductQuadrature);
    LUA_CMACRO1(GAUSS_LEGENDRE,             1);
    LUA_CMACRO1(GAUSS_CHEBYSHEV,            2);
    LUA_CMACRO1(GAUSS_LEGENDRE_LEGENDRE,    3);
    LUA_CMACRO1(GAUSS_LEGENDRE_CHEBYSHEV,   4);
    LUA_CMACRO1(CUSTOM_QUADRATURE,          5);
  LUA_FMACRO1(chiCreateLineQuadrature);
  LUA_FMACRO1(chiGetProductQuadrature);

  LUA_FMACRO1(chiCreateSLDFESQAngularQuadrature);
  LUA_FMACRO1(chiLocallyRefineSLDFESQAngularQuadrature);
  LUA_FMACRO1(chiPrintToPythonSLDFESQAngularQuadrature);
}