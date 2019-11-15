#include "ChiLua/chi_lua.h"

#include "../lbs_linear_boltzman_solver.h"
#include "ChiPhysics/chi_physics.h"
extern ChiPhysics chi_physics_handler;

//###################################################################
/**Initializes the solver.

\param SolverIndex int Handle to the solver.
 \ingroup LuaNPT
 */
int chiLBSInitialize(lua_State *L)
{
  int solver_index = lua_tonumber(L,1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzman::Solver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzman::Solver))
    {
      solver = (LinearBoltzman::Solver*)(psolver);
    }
    else
    {
      fprintf(stderr,"ERROR: Incorrect solver-type"
                     "in chiLBSInitialize\n");
      exit(EXIT_FAILURE);
    }
  }
  catch(const std::out_of_range& o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiLBSInitialize\n");
    exit(EXIT_FAILURE);
  }

  solver->Initialize();

  return 0;
}
