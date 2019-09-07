#include "../../../CHI_LUA/chi_lua.h"
#include "../../chi_math.h"
#include "../product_quadrature.h"

#include <chi_log.h>

extern ChiMath    chi_math_handler;
extern ChiLog     chi_log;

//########################################################## Create empty system
/** Creates a Product-quadrature.
 *
\param QuadratureType int Quadrature identifier.
\param Np int Number of polar angles per octant.
\param Na int Number of Azimuthal angles per octant.

##_

###QuadratureType:\n
GAUSS_LEGENDRE\n
 Gauss-Legendre quadrature for the polar angles and no quadrature rule
 for the azimuthal angle. Suitable only for 1D simulations. \n\n

GAUSS_LEGENDRE_LEGENDRE\n
 Gauss-Legendre quadrature for both the polar and azimuthal dimension.\n\n

GAUSS_LEGENDRE_CHEBYSHEV\n
 Gauss-Legendre quadrature for the polar angle but Gauss-Chebyshev
 for the azimuthal angle.\n\n

\return Returns a unique handle to the created product quadrature rule

\ingroup LuaQuadrature
\author Jan*/
int chiCreateProductQuadrature(lua_State *L)
{
  int num_args = lua_gettop(L);
  //============================================= Parse argument
  int ident = lua_tonumber(L,1);
  bool verbose = false;



  if (ident == GAUSS_LEGENDRE) //GAUSS_LEGENDRE_LEGENDRE
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
  else if (ident == GAUSS_CHEBYSHEV) //GAUSS_CHEBYSHEV
  {
    if (num_args<2)
      LuaPostArgAmountError("chiCreateProductQuadrature",2,num_args);

    int Na = lua_tonumber(L,2);
    if (num_args == 3)
      verbose = lua_toboolean(L,3);

    chi_math::ProductQuadrature* new_quad = new chi_math::ProductQuadrature;
    new_quad->InitializeWithGC(Na,verbose);
    chi_math_handler.product_quadratures.push_back(new_quad);
    int index = chi_math_handler.product_quadratures.size()-1;
    lua_pushnumber(L,index);

    if (verbose)
    {
      chi_log.Log(LOG_0)
        << "Created Gauss-Chebyshev Quadrature with "
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
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "In call to chiCreateProductQuadrature. Unsupported quadrature type"
         " supplied. Given: " << ident;
    exit(EXIT_FAILURE);
  }
  return 0;
}