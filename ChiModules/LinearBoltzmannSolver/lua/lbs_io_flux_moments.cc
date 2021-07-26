#include "ChiLua/chi_lua.h"

#include "../lbs_linear_boltzmann_solver.h"

#include "chi_log.h"
extern ChiLog&     chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Writes the flux-moments of a LBS solution to file.

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
  chi_physics::Solver* psolver;
  LinearBoltzmann::Solver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    solver = dynamic_cast<LinearBoltzmann::Solver*>(psolver);

    if (not solver)
    {
      chi_log.Log(LOG_ALLERROR)
        << __FUNCTION__ << ": Incorrect solver-type."
                           " Cannot cast to LinearBoltzmann::Solver\n";
      exit(EXIT_FAILURE);
    }
  }
  catch(const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to " << __FUNCTION__;
    exit(EXIT_FAILURE);
  }

  std::vector<double> source_moments;
  solver->MakeSourceMomentsFromPhi(source_moments);
  solver->WriteFluxMoments(file_base, source_moments);

  return 0;
}

//###################################################################
/**Reads the fluxes-moments from a file to a specific flux_moment_vector.

\param SolverIndex int Handle to the solver for which the group
is to be created.

\param file_base string Path+Filename_base to use for the output. Each location
                        will append its id to the back plus an extension ".data"

*/
int chiLBSReadFluxMoments(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if ((num_args != 2) and (num_args != 3) and (num_args != 4))
    LuaPostArgAmountError(__FUNCTION__,2,num_args);

  LuaCheckNilValue(__FUNCTION__,L,1);
  LuaCheckNilValue(__FUNCTION__,L,2);

  int      solver_index = lua_tonumber(L,1);
  std::string file_base = lua_tostring(L,2);

  bool single_file_flag = false;
  if (num_args >= 3)
  {
    LuaCheckBoolValue(__FUNCTION__, L, 3);
    single_file_flag = lua_toboolean(L, 3);
  }

  bool file_is_flux = false;
  if (num_args == 4)
  {
    LuaCheckBoolValue(__FUNCTION__, L, 4);
    file_is_flux = lua_toboolean(L, 4);
  }

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmann::Solver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    solver = dynamic_cast<LinearBoltzmann::Solver*>(psolver);

    if (not solver)
    {
      chi_log.Log(LOG_ALLERROR)
        << __FUNCTION__ << ": Incorrect solver-type."
                           " Cannot cast to LinearBoltzmann::Solver\n";
      exit(EXIT_FAILURE);
    }
  }
  catch(const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to " << __FUNCTION__;
    exit(EXIT_FAILURE);
  }

  solver->ReadFluxMoments(file_base,
                          solver->ext_src_moments_local,
                          single_file_flag);

  if (file_is_flux)
  {
    chi_log.Log() << "Making source moments from flux file." << single_file_flag << file_is_flux;
    auto temp_phi = solver->phi_old_local;
    solver->phi_old_local = solver->ext_src_moments_local;
    solver->MakeSourceMomentsFromPhi(solver->ext_src_moments_local);
    solver->phi_old_local.assign(solver->phi_old_local.size(),0.0);
  }

  return 0;
}