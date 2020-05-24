#include "../../../ChiLua/chi_lua.h"
#include "../../chi_math.h"

#include <chi_log.h>

extern ChiMath&     chi_math_handler;
extern ChiLog&     chi_log;

//########################################################## Get product quadrature
/** Get the values of a product quadrature

\param QuadHandle int Handle to an existing product quadrature.

\return Table A lua table with each entry being another table with entries
       .weight .polar .azimuthal.

\ingroup LuaQuadrature
\author Jan*/
int chiGetProductQuadrature(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiGetProductQuadrature",1,num_args);

  int handle = lua_tonumber(L,1);

  chi_math::ProductQuadrature* quad;
  try{
    quad = chi_math_handler.product_quadratures.at(handle);
  }
  catch (const std::out_of_range& o){
    chi_log.Log(LOG_ALLERROR)
      << "chiGetProductQuadrature: Invalid quadrature handle.";
    exit(EXIT_FAILURE);
  }

  lua_newtable(L);
  for (size_t n=0; n<quad->weights.size(); ++n)
  {
    lua_pushnumber(L,n+1);
    lua_newtable(L);

    lua_pushstring(L,"weight");
    lua_pushnumber(L,quad->weights[n]);
    lua_settable(L,-3);

    lua_pushstring(L,"polar");
    lua_pushnumber(L,quad->abscissae[n]->theta);
    lua_settable(L,-3);

    lua_pushstring(L,"azimuthal");
    lua_pushnumber(L,quad->abscissae[n]->phi);
    lua_settable(L,-3);

    lua_settable(L,-3);
  }

  return 1;
}