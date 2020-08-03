#include "ChiLua/chi_lua.h"

#include "../sldfe_sq.h"

#include "ChiMath/chi_math.h"
extern ChiMath&     chi_math_handler;

//###################################################################
/** Creates a Simplified Linear Discontinuous Finite Element (SLDFE)
quadrature based on Spherical Quadrilaterals (SQ). Hence SLDFE-SQ.
\param initial_refinement_level int Initial refinement level, \f$n\f$ to
       be used. The total number of angles will be \f$ 8{\times}12(n+1)^2 \f$.

\ingroup LuaQuadrature
\author Jan */
int chiCreateSLDFESQAngularQuadrature(lua_State* L)
{
  auto sldfesq = new chi_math::SimplifiedLDFESQ::Quadrature;

  int init_refinement_level = 0;

  sldfesq->GenerateInitialRefinement(init_refinement_level);


  std::shared_ptr<chi_math::AngularQuadrature> new_ang_quad =
    std::shared_ptr<chi_math::SimplifiedLDFESQ::Quadrature>(sldfesq);

  chi_math_handler.angular_quadratures.push_back(new_ang_quad);
  int index = chi_math_handler.angular_quadratures.size() - 1;
  lua_pushnumber(L,index);

  return 1;
}