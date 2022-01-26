#include"ChiLua/chi_lua.h"

#include"../Solver/diffusion_solver.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include"ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

//#############################################################################
/** Sets a property of a Diffusion solver. Please also consult the whitepaper
 * for the Diffusion solver (<a
 * href="../../whitepages/DiffusionSolver/DiffusionSolver.pdf">
 * Diffusion Whitepaper</a>)

\n\n Additional basic options can be set as indicated in \ref LuaDiffusionBasicOptions

\param SolverHandle int Handle to an existing diffusion solver.
\param PropertyName string Name for a specific property.
\param Values varying Number of inputs associated with the index.<br>

##_

###PropertyName\n
"boundary_type"\n
 Boundary type. Expects boundary index then <B>BoundaryTypeName</B>
 then type value.\n\n

\code
chiDiffusionSetProperty(solver,"boundary_type",2,"dirichlet",1.0)
\endcode

### BoundaryTypeName
reflecting\n
 Reflecting boundary conditions. Synonymous with Neumann with a
 derivative of 0.0.
             \f[ -D \hat{n}\cdot \nabla \phi = 0 \f]\n\n

dirichlet\n
 Constant value boundary condition.
 Expects to be followed by a value \f$ f \f$ associated with \f$ \phi \f$.
            \f[ \phi = f \f]\n\n

neumann\n
 Constant derivative boundary condition. Expects to be followed
 by a constant \f$ f \f$ representing
                    \f[ -D \hat{n}\cdot \nabla \phi = f \f]\n\n

vacuum\n
 Vacuum boundary conditions. More appropriate to neutron diffusion.
   \f[ \frac{1}{4}\phi + \frac{1}{2} D \hat{n}\cdot \nabla \phi = 0 \f]\n\n

robin\n
 Robin boundary condition of the form
                   \f[ a \phi + b D \hat{n}\cdot \nabla \phi = f \f]\n\n

\ingroup LuaDiffusion
\author Jan*/
int chiDiffusionSetProperty(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 2)
    LuaPostArgAmountError(fname, num_args, 2);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  //============================================= Get solver
  LuaCheckNumberValue(fname, L, 1);
  int solver_index = lua_tonumber(L,1);
  chi_diffusion::Solver* solver;

  try{
    solver = (chi_diffusion::Solver*)chi_physics_handler.solver_stack.at(solver_index);
  }
  catch(const std::out_of_range& o){
    chi_log.Log(LOG_0ERROR) << "ERROR: Invalid solver handle." << std::endl;
    exit(EXIT_FAILURE);
  }

  //============================================= Get property index
  LuaCheckStringValue(fname, L, 2);
  const std::string property_name = lua_tostring(L, 2);

  //============================================= Handle properties
  if (property_name == "boundary_type")
  {
    if (num_args < 4)
    {
      chi_log.Log(LOG_0ERROR)
      << "Invalid amount of arguments used in"
      << " chiDiffusionSetProperty(...,\"boundary_type\".... "
      << " At least 4 arguments are expected.";
      exit(EXIT_FAILURE);
    }
    LuaCheckNumberValue(fname, L, 3);
    const int bound_index = lua_tonumber(L,3);

    LuaCheckStringValue(fname, L, 4);
    const std::string type_name = lua_tostring(L, 4);

    if (type_name == "reflecting")
    {
      if (num_args != 4)
      {
        chi_log.Log(LOG_0ERROR)
          << "Invalid amount of arguments used in"
          << " chiDiffusionSetProperty(...,\"boundary_type\","
          << bound_index << ",\"reflecting\". "
          << " 4 arguments are expected.";
        exit(EXIT_FAILURE);
      }

      chi_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = chi_diffusion::BoundaryType::Reflecting;

      solver->boundary_preferences.insert(std::make_pair(bound_index,bndry_info));

      chi_log.Log(LOG_0) << "Boundary " << bound_index << " set as "
                         << "Reflecting.";
    }
    else if (type_name == "dirichlet")
    {
      if (num_args != 5)
      {
        chi_log.Log(LOG_0ERROR)
          << "Invalid amount of arguments used in"
          << " chiDiffusionSetProperty(...,\"boundary_type\","
          << bound_index << ",\"dirichlet\". "
          << " 5 arguments are expected.";
        exit(EXIT_FAILURE);
      }
      LuaCheckNumberValue(fname, L, 5);
      double b_value = lua_tonumber(L,5);

      chi_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = chi_diffusion::BoundaryType::Dirichlet;
      bndry_info.second = {b_value};
      solver->boundary_preferences.insert(std::make_pair(bound_index,bndry_info));

      chi_log.Log(LOG_0) << "Boundary " << bound_index << " set as "
                         << "Dirichlet with value " << b_value;
    }
    else if (type_name == "neumann")
    {
      if (num_args != 5)
      {
        chi_log.Log(LOG_0ERROR)
          << "Invalid amount of arguments used in"
          << " chiDiffusionSetProperty(...,\"boundary_type\","
          << bound_index << ",\"neumann\". "
          << " 5 arguments are expected.";
        exit(EXIT_FAILURE);
      }
      LuaCheckNumberValue(fname, L, 5);
      double f_value = lua_tonumber(L,5);

      chi_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = chi_diffusion::BoundaryType::Robin;
      bndry_info.second = {0.0,1.0,f_value};
      solver->boundary_preferences.insert(std::make_pair(bound_index,bndry_info));

      chi_log.Log(LOG_0) << "Boundary " << bound_index << " set as "
                         << "Neumann with f = ("
                         << f_value << ") ";
    }
    else if (type_name == "vacuum")
    {
      if (num_args != 4)
      {
        chi_log.Log(LOG_0ERROR)
          << "Invalid amount of arguments used in"
          << " chiDiffusionSetProperty(...,\"boundary_type\","
          << bound_index << ",\"vacuum\". "
          << " 4 arguments are expected.";
        exit(EXIT_FAILURE);
      }

      chi_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = chi_diffusion::BoundaryType::Robin;
      bndry_info.second = {0.25,0.5,0.0};
      solver->boundary_preferences.insert(std::make_pair(bound_index,bndry_info));

      chi_log.Log(LOG_0) << "Boundary " << bound_index << " set as "
                         << "Vacuum.";
    }
    else if (type_name == "robin")
    {
      if (num_args != 7)
      {
        chi_log.Log(LOG_0ERROR)
          << "Invalid amount of arguments used in"
          << " chiDiffusionSetProperty(...,\"boundary_type\","
          << bound_index << ",\"robin\". "
          << " 7 arguments are expected.";
        exit(EXIT_FAILURE);
      }
      LuaCheckNumberValue(fname, L, 5);
      LuaCheckNumberValue(fname, L, 6);
      LuaCheckNumberValue(fname, L, 7);

      double a_value = lua_tonumber(L,5);
      double b_value = lua_tonumber(L,6);
      double f_value = lua_tonumber(L,7);

      chi_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = chi_diffusion::BoundaryType::Robin;
      bndry_info.second = {a_value,b_value,f_value};
      solver->boundary_preferences.insert(std::make_pair(bound_index,bndry_info));

      chi_log.Log(LOG_0) << "Boundary " << bound_index << " set as "
                         << "Robin with a,b,f = ("
                         << a_value << ","
                         << b_value << ","
                         << f_value << ") ";
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Unsupported boundary type encountered in call to "
        << "chiDiffusionSetProperty(..,\"boundary_type\",.. :"
        << type_name;
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    chi_log.Log(LOG_0ERROR) << "Invalid property in chiDiffusionSetProperty.";
    exit(EXIT_FAILURE);
  }
  return 0;
}
