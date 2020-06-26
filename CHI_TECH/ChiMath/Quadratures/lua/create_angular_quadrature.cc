#include "../../../ChiLua/chi_lua.h"
#include "ChiMath/chi_math.h"
#include "ChiMath/Quadratures/angular_quadrature_base.h"

extern ChiMath&     chi_math_handler;

#include <chi_log.h>

extern ChiLog& chi_log;

//########################################################## Create empty system
/** Creates an angular quadrature.
 *
\param azimuthal_angles array A lua table with N entries each being an azimuthal
                              angle.
\param polar_angles array A lua table with N entries each being a polar angle
                              angle.
\param weight array A lua table with N entries each being a quadrature weight.


\return Returns a unique handle to the created angular quadrature

\ingroup LuaQuadrature
\author Jan*/
int chiCreateCustomAngularQuadrature(lua_State *L)
{
  size_t num_args = lua_gettop(L);

  if (num_args != 3)
    LuaPostArgAmountError("chiCreateCustomAngularQuadrature",3,num_args);

  LuaCheckNilValue("chiCreateCustomAngularQuadrature",L,1);
  LuaCheckNilValue("chiCreateCustomAngularQuadrature",L,2);
  LuaCheckNilValue("chiCreateCustomAngularQuadrature",L,3);

  LuaCheckTableValue("chiCreateCustomAngularQuadrature",L,1);
  LuaCheckTableValue("chiCreateCustomAngularQuadrature",L,2);
  LuaCheckTableValue("chiCreateCustomAngularQuadrature",L,3);

  int Na = lua_rawlen(L,1);
  int Np = lua_rawlen(L,2);
  int Nw = lua_rawlen(L,3);

  if ((Na-Np != 0) or (Na-Nw !=0))
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiCreateCustomAngularQuadrature: Tables lengths supplied "
         "are not of equal lengths.";
    exit(EXIT_FAILURE);
  }

  std::vector<double> azi_angles(Na,0.0);
  std::vector<double> pol_angles(Na,0.0);
  std::vector<double> weights(Na,0.0);

  for (int n=1; n<=Na; ++n)
  {
    lua_pushnumber(L,n);
    lua_gettable(L,1);
    azi_angles[n-1] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  for (int n=1; n<=Na; ++n)
  {
    lua_pushnumber(L,n);
    lua_gettable(L,2);
    pol_angles[n-1] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  for (int n=1; n<=Na; ++n)
  {
    lua_pushnumber(L,n);
    lua_gettable(L,3);
    weights[n-1] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  chi_log.Log(LOG_0) << "Creating Custom Angular Quadrature\n";

  auto angular_quadrature = std::make_shared<chi_math::AngularQuadrature>();
  angular_quadrature->InitializeWithCustom(azi_angles,pol_angles,weights);

  chi_math_handler.angular_quadratures.push_back(angular_quadrature);
  int index = chi_math_handler.angular_quadratures.size()-1;
  lua_pushnumber(L,index);

  return 1;
}