#include "../../ChiLua/chi_lua.h"
#include<iostream>
#include "../chi_physics.h"
#include "../../ChiMesh/Region/chi_region.h"
#include "../../ChiMesh/MeshHandler/chi_meshhandler.h"

extern ChiPhysics chi_physics_handler;

/** \defgroup LuaSolver Solvers
 * \ingroup LuaPhysics*/

//#############################################################################
/** Adds a region to a solver.
 *
\param SolverHandle int Handle to the solver.
\param RegionHandle int Handle to the region.

\ingroup LuaSolver
\author Jan*/
int chiSolverAddRegion(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  int solver_index = lua_tonumber(L,1);
  int region_index = lua_tonumber(L,2);

  //======================================================= Getting solver
  chi_physics::Solver* solver;
  try{
    solver = chi_physics_handler.solver_stack.at(solver_index);
  }
  catch(std::out_of_range o)
  {
    std::cout << "Invalid solver handle" << std::endl;
    return 0;
  }

  //======================================================= Getting region
  chi_mesh::Region* region;
  try{
    region = cur_hndlr->region_stack.at(region_index);
  }
  catch(std::out_of_range o)
  {
    std::cout << "Invalid region handle" << std::endl;
    return 0;
  }

  solver->AddRegion(region);

  return 0;
}
