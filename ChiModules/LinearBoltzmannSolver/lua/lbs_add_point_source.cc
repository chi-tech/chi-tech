#include "chi_lua.h"
#include "lbs_lua_utils.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Adds a point source to an LBS solver.
\param SolverIndex int Handle to the solver.
\param Location_x double X-location.
\param Location_y double Y-location.
\param Location_z double Z-location.
\param Strength table Source strength as a multigroup vector.

 \ingroup LuaLBS
 */
int chiLBSAddPointSource(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 5)
    LuaPostArgAmountError(fname, 5, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  LuaCheckNilValue(fname, L, 4);
  LuaCheckNilValue(fname, L, 5);

  //============================================= Get pointer to solver
  const int solver_index = lua_tonumber(L,1);
  auto lbs_solver = lbs::lua_utils::GetSolverByHandle(solver_index, fname);

  //============================================= Get other arguments
  const double x = lua_tonumber(L, 2);
  const double y = lua_tonumber(L, 3);
  const double z = lua_tonumber(L, 4);

  const chi_mesh::Vector3 location(x,y,z);


  LuaCheckTableValue(fname, L, 5);

  std::vector<double> strength;
  LuaPopulateVectorFrom1DArray(fname,L, 5, strength);

  lbs_solver->point_sources.emplace_back(location, strength);

  chi_log.Log() << "LBS: Added point source at "
                << location.PrintStr();

  return 0;
}