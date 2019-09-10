#include <ChiLua/chi_lua.h>

#include"../Solver/solver_montecarlon.h"

#include <ChiPhysics/chi_physics.h>

extern ChiPhysics chi_physics_handler;

#include <chi_log.h>

//#############################################################################
/** Executes a MonteCarlon solver.

\param SolverHandle int Handle to the montecarlo solver.
\param PropertyIndex int Code for a specific property.



##_

###PropertyIndex\n
MC_NUM_PARTICLES\n
 Number of particles to run. Expects to be followed by an integer specifying
 the amount of particles to run. Default 1000.\n\n

MC_TFC_UPDATE_DIV\n
 Number of divisions of the number of particles to use for tally fluctuation
 chart (TFC) update. Expects to be followed by an integer specifying the number
 of bins. Default 10.\n\n

MC_MONOENERGETIC\n
 Forces the scattering out of a group to be treated like absorbtion.
 Expects to be followed by a boolean value. Default false.\n\n

MC_SCATTERING_ORDER\n
 Sets the scattering order used for building discrete scattering angles.
 Expect to be followed by an integer specifying the order. Default 10.
 Note: when this number is set greater than the scattering order available
 in the provided cross-sections then the scattering order will default to that
 available.\n\n

MC_FORCE_ISOTROPIC\n
 Flag forcing isotropic scattering. Expects to be followed by a boolean value.
 Default false.\n\n

MC_TALLY_MULTIPLICATION_FACTOR\n
 Classical global tally multiplication factor to be applied after normalization
 per source particle. Expects to be followed by a float. Default 1.0.\n\n

\ingroup LuaMonteCarlon
\author Jan*/
int chiMonteCarlonSetProperty(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args < 3)
    LuaPostArgAmountError("chiMonteCarlonSetProperty",3,num_args);

  chi_physics::Solver* solver = nullptr;
  try{
    solver = chi_physics_handler.solver_stack.at(lua_tonumber(L,1));
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiMonteCarlonSetProperty: Invalid solver handle. "
      << lua_tonumber(L,1);
    exit(EXIT_FAILURE);
  }

  chi_montecarlon::Solver* mcsolver;
  if (typeid(*solver) == typeid(chi_montecarlon::Solver))
    mcsolver = (chi_montecarlon::Solver*)solver;
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiMonteCarlonSetProperty: Solver pointed to by solver handle is "
      << " not a MonteCarlo solver.";
    exit(EXIT_FAILURE);
  }

  //============================================= Process property index
  int property_index = lua_tonumber(L,2);
  if (property_index == MC_NUM_PARTICLES)
  {
    unsigned long long num_part = lua_tonumber(L,3);

    mcsolver->num_particles = num_part;
  }
  else if (property_index == MC_TFC_UPDATE_INTVL)
  {
    int tfc_updt_intvl = lua_tonumber(L,3);

    mcsolver->tfc_update_interval = tfc_updt_intvl;
  }
  else if (property_index == MC_SCATTERING_ORDER)
  {
    int scatorder = lua_tonumber(L,3);

    mcsolver->scattering_order = scatorder;
  }
  else if (property_index == MC_MONOENERGETIC)
  {
    bool mono = lua_toboolean(L,3);

    mcsolver->mono_energy = mono;
  }
  else if (property_index == MC_FORCE_ISOTROPIC)
  {
    bool iso = lua_toboolean(L,3);

    mcsolver->force_isotropic = iso;
  }
  else if (property_index == MC_GROUP_BOUNDS)
  {
    int hi = lua_tonumber(L,3);
    int lo = lua_tonumber(L,3);

    mcsolver->group_hi_bound = hi;
    mcsolver->group_lo_bound = lo;
  }
  else if (property_index == MC_TALLY_MERGE_INTVL)
  {
    unsigned long long tal_merg_invtl = lua_tonumber(L,3);

    mcsolver->tally_rendezvous_intvl = tal_merg_invtl;
  }
  else if (property_index == MC_TALLY_MULTIPLICATION_FACTOR)
  {
    double tmf = lua_tonumber(L,3);

    mcsolver->tally_multipl_factor = tmf;
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiMonteCarlonSetProperty:Invalid property index supplied. "
      << property_index;
    exit(EXIT_FAILURE);
  }

  return 0;
}