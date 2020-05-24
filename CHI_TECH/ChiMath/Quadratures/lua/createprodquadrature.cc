#include "../../../ChiLua/chi_lua.h"
#include "../../chi_math.h"

#include <chi_log.h>

extern ChiMath&     chi_math_handler;
extern ChiLog&     chi_log;

//########################################################## Create empty system
/** Creates a Product-quadrature.
 *
\param QuadratureType int Quadrature identifier.
\param values varying Varying options based on the quadrature type.

##_

###QuadratureType:\n
GAUSS_LEGENDRE\n
 Gauss-Legendre quadrature for the polar angles and no quadrature rule
 for the azimuthal angle. Suitable only for 1D simulations. Expects
 to be followed by the number of angles Np. Optionally a verbosity flag
 can be added.\n\n

GAUSS_LEGENDRE_LEGENDRE\n
 Gauss-Legendre quadrature for both the polar and azimuthal dimension.
 Expects to be followed by number of Azimuthal and Polar angles.
 Optionally a verbosity flag can be added.\n\n

GAUSS_LEGENDRE_CHEBYSHEV\n
 Gauss-Legendre quadrature for the polar angle but Gauss-Chebyshev
 for the azimuthal angle.
 Expects to be followed by number of Azimuthal and Polar angles.
 Optionally a verbosity flag can be added.\n\n

CUSTOM_QUADRATURE\n
 Expects to be followed by three lua tables. The first table is an array,
 of length Na, of the azimuthal angles (radians). The second table is an array,
 of length Np, of the polar angles (radians). The third table is an array, of
 length Na*Np, and contains the weight associated with each angle pair.
 Optionally a verbosity flag can be added.\n\n


\return Returns a unique handle to the created product quadrature rule

\ingroup LuaQuadrature
\author Jan*/
int chiCreateProductQuadrature(lua_State *L)
{
  int num_args = lua_gettop(L);
  //============================================= Parse argument
  int ident = lua_tonumber(L,1);
  bool verbose = false;



  if (ident == GAUSS_LEGENDRE) //GAUSS_LEGENDRE
  {
    if (num_args<2)
      LuaPostArgAmountError("chiCreateProductQuadrature",2,num_args);

    int Np = lua_tonumber(L,2);
    if (num_args == 3)
      verbose = lua_toboolean(L,3);

    chi_math::ProductQuadrature* new_quad = new chi_math::ProductQuadrature;
    new_quad->InitializeWithGL(Np,verbose);
    chi_math_handler.product_quadratures.push_back(new_quad);
    int index = chi_math_handler.product_quadratures.size()-1;
    lua_pushnumber(L,index);

    if (verbose)
    {
      chi_log.Log(LOG_0)
        << "Created Gauss-Legendre Quadrature with "
        << new_quad->azimu_ang.size()
        << " azimuthal angles and "
        << new_quad->polar_ang.size()
        << " polar angles.";
    }

    return 1;
  }
  else if (ident == GAUSS_LEGENDRE_LEGENDRE) //GAUSS_LEGENDRE_LEGENDRE
  {
    if (num_args<3)
      LuaPostArgAmountError("chiCreateProductQuadrature",3,num_args);

    int Np = lua_tonumber(L,2);
    int Na = lua_tonumber(L,3);
    if (num_args == 4)
      verbose = lua_toboolean(L,4);

    chi_math::ProductQuadrature* new_quad = new chi_math::ProductQuadrature;
    new_quad->InitializeWithGLL(Np,Na,verbose);
    chi_math_handler.product_quadratures.push_back(new_quad);
    int index = chi_math_handler.product_quadratures.size()-1;
    lua_pushnumber(L,index);

    if (verbose)
    {
      chi_log.Log(LOG_0)
        << "Created Gauss-Legendre-Legendre Quadrature with "
        << new_quad->azimu_ang.size()
        << " azimuthal angles and "
        << new_quad->polar_ang.size()
        << " polar angles.";
    }

    return 1;
  }
  else if (ident == GAUSS_LEGENDRE_CHEBYSHEV) //GAUSS_LEGENDRE_CHEBYSHEV
  {
    if (num_args<3)
      LuaPostArgAmountError("chiCreateProductQuadrature",3,num_args);

    int Np = lua_tonumber(L,2);
    int Na = lua_tonumber(L,3);
    if (num_args == 4)
      verbose = lua_toboolean(L,4);

    chi_math::ProductQuadrature* new_quad = new chi_math::ProductQuadrature;
    new_quad->InitializeWithGLC(Np,Na,verbose);
    chi_math_handler.product_quadratures.push_back(new_quad);
    int index = chi_math_handler.product_quadratures.size()-1;
    lua_pushnumber(L,index);

    if (verbose)
    {
      chi_log.Log(LOG_0)
      << "Created Gauss-Legendre-Chebyshev Quadrature with "
      << new_quad->azimu_ang.size()
      << " azimuthal angles and "
      << new_quad->polar_ang.size()
      << " polar angles.";
    }

    return 1;
  }
  else if (ident == CUSTOM_QUADRATURE) //CUSTOM_QUADRATURE
  {
    if (num_args<4)
      LuaPostArgAmountError("chiCreateProductQuadrature:CUSTOM_QUADRATURE",3,num_args);

    if (not lua_istable(L,2))
    {
      chi_log.Log(LOG_ALLERROR)
        << "chiCreateProductQuadrature:CUSTOM_QUADRATURE, second argument must "
        << "be a lua table.";
      exit(EXIT_FAILURE);
    }
    if (not lua_istable(L,3))
    {
      chi_log.Log(LOG_ALLERROR)
        << "chiCreateProductQuadrature:CUSTOM_QUADRATURE, third argument must "
        << "be a lua table.";
      exit(EXIT_FAILURE);
    }
    if (not lua_istable(L,4))
    {
      chi_log.Log(LOG_ALLERROR)
        << "chiCreateProductQuadrature:CUSTOM_QUADRATURE, fourth argument must "
        << "be a lua table.";
      exit(EXIT_FAILURE);
    }
    if (num_args == 5)
      verbose = lua_toboolean(L,4);

    int Na = lua_rawlen(L,2);
    int Np = lua_rawlen(L,3);
    int Nw = lua_rawlen(L,4);

    std::vector<double> azimuthal(Na,0.0);
    std::vector<double> polar(Np,0.0);
    std::vector<double> weights(Nw,0.0);

    for (int n=1; n<=Na; ++n)
    {
      lua_pushnumber(L,n);
      lua_gettable(L,2);
      azimuthal[n-1] = lua_tonumber(L,-1);
      lua_pop(L,1);
    }
    for (int n=1; n<=Np; ++n)
    {
      lua_pushnumber(L,n);
      lua_gettable(L,3);
      polar[n-1] = lua_tonumber(L,-1);
      lua_pop(L,1);
    }
    for (int n=1; n<=Nw; ++n)
    {
      lua_pushnumber(L,n);
      lua_gettable(L,4);
      weights[n-1] = lua_tonumber(L,-1);
      lua_pop(L,1);
    }

    chi_log.Log(LOG_0) << Na << " " << Np << " " << Nw;

    auto new_quad = new chi_math::ProductQuadrature;
    new_quad->InitializeWithCustom(azimuthal,polar,weights,verbose);
    chi_math_handler.product_quadratures.push_back(new_quad);
    int index = chi_math_handler.product_quadratures.size()-1;
    lua_pushnumber(L,index);

    if (verbose)
    {
      chi_log.Log(LOG_0)
        << "Created Custom Quadrature with "
        << new_quad->azimu_ang.size()
        << " azimuthal angles and "
        << new_quad->polar_ang.size()
        << " polar angles.";
    }

    return 1;
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "In call to chiCreateProductQuadrature. Unsupported quadrature type"
         " supplied. Given: " << ident;
    exit(EXIT_FAILURE);
  }
  return 0;
}