#include "ChiLua/chi_lua.h"

#include "../lbs_linear_boltzmann_solver.h"

#include "lbs_lua_utils.h"

#include "chi_log.h"
extern ChiLog&     chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Writes the flux-moments of a LBS solution to file (phi_old_local).

\param SolverIndex int Handle to the solver for which the group
is to be created.

\param file_base string Path+Filename_base to use for the output. Each location
                        will append its id to the back plus an extension ".data"

*/
int chiLBSWriteFluxMoments(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(__FUNCTION__,2,num_args);

  LuaCheckNilValue(__FUNCTION__,L,1);
  LuaCheckNilValue(__FUNCTION__,L,2);

  int      solver_index = lua_tonumber(L,1);
  std::string file_base = lua_tostring(L,2);

  //============================================= Get pointer to solver
  auto solver =
    LinearBoltzmann::lua_utils::GetSolverByHandle(solver_index, __FUNCTION__);

  solver->WriteFluxMoments(file_base, solver->phi_old_local);

  return 0;
}

//###################################################################
/**Creates scattered source-moments, based on a LBS solution, and writes them
 * to file.

\param SolverIndex int Handle to the solver for which the group
is to be created.

\param file_base string Path+Filename_base to use for the output. Each location
                        will append its id to the back plus an extension ".data"

*/
int chiLBSCreateAndWriteSourceMoments(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(__FUNCTION__,2,num_args);

  LuaCheckNilValue(__FUNCTION__,L,1);
  LuaCheckNilValue(__FUNCTION__,L,2);

  int      solver_index = lua_tonumber(L,1);
  std::string file_base = lua_tostring(L,2);

  //============================================= Get pointer to solver
  auto solver =
    LinearBoltzmann::lua_utils::GetSolverByHandle(solver_index, __FUNCTION__);

  auto source_moments = solver->MakeSourceMomentsFromPhi();
  solver->WriteFluxMoments(file_base, source_moments);

  return 0;
}

//###################################################################
/**Reads flux-moments from a file and creates a scattering source from these
 * moments to be used instead of a regular material/boundary source.

\param SolverIndex int Handle to the solver for which the group
is to be created.

\param file_base string Path+Filename_base to use for the output. Each location
                        will append its id to the back plus an extension ".data"

\param single_file_flag bool (Optional) Flag indicating that the file is a
                             single stand-alone file. The file_base will then
                             be used without adding the location-id, but still
                             with the ".data" appended. Default: false.

*/
int chiLBSReadFluxMomentsAndMakeSourceMoments(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if ((num_args != 2) and (num_args != 3))
    LuaPostArgAmountError(__FUNCTION__,2,num_args);

  LuaCheckNilValue(__FUNCTION__,L,1);
  LuaCheckNilValue(__FUNCTION__,L,2);

  int      solver_index = lua_tonumber(L,1);
  std::string file_base = lua_tostring(L,2);

  bool single_file_flag = false;
  if (num_args == 3)
  {
    LuaCheckBoolValue(__FUNCTION__, L, 3);
    single_file_flag = lua_toboolean(L, 3);
  }

  //============================================= Get pointer to solver
 auto solver =
   LinearBoltzmann::lua_utils::GetSolverByHandle(solver_index, __FUNCTION__);

  solver->ReadFluxMoments(file_base,
                          solver->ext_src_moments_local,
                          single_file_flag);

  chi_log.Log() << "Making source moments from flux file.";
  auto temp_phi = solver->phi_old_local;
  solver->phi_old_local = solver->ext_src_moments_local;
  solver->ext_src_moments_local = solver->MakeSourceMomentsFromPhi();
  solver->phi_old_local = temp_phi;

  return 0;
}

//###################################################################
/**Reads the source-moments from a file to a specific ext_src_moments_local-vector
 * to be used instead of a regular material/boundary source.

\param SolverIndex int Handle to the solver for which the group
is to be created.

\param file_base string Path+Filename_base to use for the output. Each location
                        will append its id to the back plus an extension ".data"

\param single_file_flag bool (Optional) Flag indicating that the file is a
                             single stand-alone file. The file_base will then
                             be used without adding the location-id, but still
                             with the ".data" appended. Default: false.
*/
int chiLBSReadSourceMoments(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if ((num_args != 2) and (num_args != 3))
    LuaPostArgAmountError(__FUNCTION__,2,num_args);

  LuaCheckNilValue(__FUNCTION__,L,1);
  LuaCheckNilValue(__FUNCTION__,L,2);

  int      solver_index = lua_tonumber(L,1);
  std::string file_base = lua_tostring(L,2);

  bool single_file_flag = false;
  if (num_args == 3)
  {
    LuaCheckBoolValue(__FUNCTION__, L, 3);
    single_file_flag = lua_toboolean(L, 3);
  }

  //============================================= Get pointer to solver
  auto solver =
    LinearBoltzmann::lua_utils::GetSolverByHandle(solver_index, __FUNCTION__);

  solver->ReadFluxMoments(file_base,
                          solver->ext_src_moments_local,
                          single_file_flag);

  return 0;
}

//###################################################################
/**Reads flux-moments from a file to phi_old_local (the initial flux solution).

\param SolverIndex int Handle to the solver for which the group
is to be created.

\param file_base string Path+Filename_base to use for the output. Each location
                        will append its id to the back plus an extension ".data"

\param single_file_flag bool (Optional) Flag indicating that the file is a
                             single stand-alone file. The file_base will then
                             be used without adding the location-id, but still
                             with the ".data" appended. Default: false.
*/
int chiLBSReadFluxMoments(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if ((num_args != 2) and (num_args != 3))
    LuaPostArgAmountError(__FUNCTION__,2,num_args);

  LuaCheckNilValue(__FUNCTION__,L,1);
  LuaCheckNilValue(__FUNCTION__,L,2);

  int      solver_index = lua_tonumber(L,1);
  std::string file_base = lua_tostring(L,2);

  bool single_file_flag = false;
  if (num_args == 3)
  {
    LuaCheckBoolValue(__FUNCTION__, L, 3);
    single_file_flag = lua_toboolean(L, 3);
  }

  //============================================= Get pointer to solver
  auto solver =
    LinearBoltzmann::lua_utils::GetSolverByHandle(solver_index, __FUNCTION__);

  solver->ReadFluxMoments(file_base,
                          solver->phi_old_local,
                          single_file_flag);

  return 0;
}