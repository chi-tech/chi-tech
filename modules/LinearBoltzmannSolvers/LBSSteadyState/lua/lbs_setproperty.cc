#include "ChiLua/chi_lua.h"
#include "lbs_lua_utils.h"

#define DISCRETIZATION_METHOD 1
  #define PWLD   3

#define BOUNDARY_CONDITION 3
  #define XMAX 31
  #define XMIN 32
  #define YMAX 33
  #define YMIN 34
  #define ZMAX 35
  #define ZMIN 36

#define SCATTERING_ORDER 4

#define SWEEP_EAGER_LIMIT 5

#define READ_RESTART_DATA 6

#define WRITE_RESTART_DATA 7

#define SAVE_ANGULAR_FLUX 8

#define USE_SOURCE_MOMENTS 9

#define VERBOSE_INNER_ITERATIONS 10

#define VERBOSE_OUTER_ITERATIONS 11

#define USE_PRECURSORS 12

#include "chi_runtime.h"
#include "chi_log.h"


//###################################################################
/**Set LBS property.
\param SolverIndex int Handle to the solver for which the set is to be created.
\param PropertyIndex int Code for a specific property.


##_

###PropertyIndex
DISCRETIZATION_METHOD\n
 Discretization method.\n\n

BOUNDARY_CONDITION\n
 Boundary condition type. See BoundaryIdentify.\n\n

SCATTERING_ORDER\n
 Defines the level of harmonic expansion for the scattering source.Default 1.
 Expects to be followed by an integer.\n\n

SWEEP_EAGER_LIMIT\n
 The eager limit to be used in message size during sweep initialization.
 This expects to be followed by a size in bytes (Max 64,0000).Default 32,000.
 See note below.\n\n

READ_RESTART_DATA\n
 Indicates the reading of restart data from restart file.
 The value can be followed by two
 optional strings. The first is the folder name which can be relative or
 absolute, and the second is the file base name. These are defaulted to
 "YRestart" and "restart" respectively.\n\n

USE_SOURCE_MOMENTS\n
 Flag for using a vector of source moments instead the regular material/boundary
  source. Default false. This expects
 to be followed by a boolean.\n\n

VERBOSE_INNER_ITERATIONS\n
 Flag for printing inner iteration information. This is primarily used
 for printing information related to group-set-level iterative methods.
 Default true. Expects to be followed by a boolean.\n\n

VERBOSE_OUTER_ITERATIONS\n
 Flag for printing outer iteration information. This is primarily used
 for printing information aggregated over group sets such as k-eigenvalue
 iterations. Default true. Expects to be followed by a boolean.\n\n

USE_PRECURSORS\n
 Flag for using delayed neutron precursors. Default false. This expects
 to be followed by a boolean.\n\n

\code
chiLBSSetProperty(phys1,READ_RESTART_DATA,"YRestart1")
\endcode

WRITE_RESTART_DATA\n
 Indicates the writing of restart data to restart files.
 The value can be followed by two optional strings and a number
 optional strings. The first string is the folder name which can be relative or
 absolute, and the second string is the file base name. The number is the time
 interval (in minutes) for a restart write to be triggered (apart from GMRES
 restarts and the conclusion of groupset completions) .These are defaulted to
 "YRestart", "restart" and 30 minutes respectively.\n\n

\code
chiLBSSetProperty(phys1,WRITE_RESTART_DATA,"YRestart1","restart",1)
\endcode

###Discretization methods
 PWLD = Piecewise Linear Finite Element.\n

###BoundaryIdentify
This value follows the argument BOUNDARY_CONDITION and identifies which
boundary is under consideration. Right now only boundaries aligned with
cartesian axes are considered. Followed by LBSBoundaryType.\n
XMAX = Right boundary \n
XMIN = Left boundary \n
YMAX = Front boundary \n
YMIN = Back boundary \n
ZMAX = Top boundary \n
ZMIN = Bottom boundary \n

###LBSBoundaryType
Specifies the type of boundary. Depending on the type this argument needs
to be followed by one or more values. Note: By default all boundaries are
type VACUUM.\n
\n
LBSBoundaryTypes.VACUUM\n
Specifies a vaccuum boundary condition. It is not followed by any value.\n
\n
\n
LBSBoundaryTypes.INCIDENT_ISOTROPIC\n
Incident isotropic flux. This argument needs to be followed by a lua table
index 1 to G where G is the amount of energy groups. Note internally this
is mapped as 0 to G-1.\n
\n
LBSBoundaryTypes.REFLECTING\n
Reflecting boundary condition. Beware, when opposing reflecting boundary
conditions are used this enduces a cyclic dependency which will increase the
iteration convergence behavior.



###Note on the Eager limit
The eager limit is the message size limit before which non-blocking MPI send
calls will execute without waiting for a matching receive call. The limit is
platform dependent but in general 64 kb. Some systems have 32 kb as a limit
and therefore we use that as a default limit in ChiTech. There is a fine
interplay between message size and the shear amount of messages that will be
sent. In general smaller messages tend to be more efficient, however, when
there are too many small messages being sent around the communication system
on the given platform will start to suffer. One can gain a small amount of
parallel efficiency by lowering this limit, however, there is a point where
the parallel efficiency will actually get worse so use with caution.

\ingroup LuaLBS*/
int chiLBSSetProperty(lua_State *L)
{
  const std::string fname = "chiLBSSetProperty";
  const int numArgs = lua_gettop(L);
  if (numArgs < 2)
    LuaPostArgAmountError(fname, 2, numArgs);

  LuaCheckNilValue(fname, L, 1);

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);
  auto& lbs_solver = chi::GetStackItem<lbs::SteadyStateSolver>(chi::solver_stack,
                                                               solver_handle,
                                                               fname);

  //============================================= Get property index
  LuaCheckNilValue(fname, L, 2);

  const int property = lua_tonumber(L,2);


  //============================================= Handle properties
  if (property == DISCRETIZATION_METHOD)
  {
    LuaCheckNilValue(fname, L, 3);

    const int method = lua_tonumber(L,3);

    typedef chi_math::SpatialDiscretizationType SDMType;

    if (method == PWLD)
      lbs_solver.options.sd_type = SDMType::PIECEWISE_LINEAR_DISCONTINUOUS;
    else
      throw std::invalid_argument(
        "Invalid option for Discretization method in chiLBSSetProperty.\n");
  }
  else if (property == BOUNDARY_CONDITION)
  {
    if (numArgs<4)
      LuaPostArgAmountError("chiLBSSetProperty",4,numArgs);

    LuaCheckNilValue(fname, L, 3);
    LuaCheckNilValue(fname, L, 4);

    const int bident = lua_tonumber(L,3);
    const int btype  = lua_tonumber(L,4);

    if (!((bident>=XMAX) && (bident<=ZMIN)))
    {
      chi::log.LogAllError()
        << "Unknown boundary identifier encountered "
           "in call to chiLBSSetProperty";
      chi::Exit(EXIT_FAILURE);
    }

    const int bid = bident - 31;

    if (btype == (int)lbs::BoundaryType::VACUUM)
    {
      lbs_solver.boundary_types[bid].first = lbs::BoundaryType::VACUUM;
      chi::log.Log() << "Boundary " << bid << " set to Vacuum.";
    }
    else if (btype == (int)lbs::BoundaryType::INCIDENT_ISOTROPIC)
    {
      if (numArgs!=5)
        LuaPostArgAmountError("chiLBSSetProperty",5,numArgs);

      if (lbs_solver.groups.empty())
      {
        chi::log.Log0Error()
          << "In call to chiLBSSetProperty, setting "
          << "incident isotropic flux boundary type: Number of solver groups"
          << " is zero. Boundary fluxes can only be set after group structure"
          << " has been defined.";
        chi::Exit(EXIT_FAILURE);
      }

      if (!lua_istable(L,5))
      {
        chi::log.LogAllError()
          << "In call to chiLBSSetProperty, setting "
          << "incident isotropic flux boundary type,"
          << " argument 5 should be a lua table and was detected as"
             " not being one.";
        chi::Exit(EXIT_FAILURE);
      }

      const size_t table_len = lua_rawlen(L,5);
      std::vector<double> values(table_len,0.0);
      for (int g=0; g<table_len; g++)
      {
        lua_pushnumber(L,g+1);
        lua_gettable(L,5);
        values[g] = lua_tonumber(L,-1);
        lua_pop(L,1);
      }

      if (table_len != lbs_solver.groups.size())
      {
        chi::log.Log0Error()
          << "In call to chiLBSSetProperty, setting "
          << "incident isotropic flux boundary type: "
          << "Number of groups in boundary flux specification is "
          << table_len << " but solver has a total of "
          << lbs_solver.groups.size() << " groups. These two must be equal.";
        chi::Exit(EXIT_FAILURE);
      }

      lbs_solver.incident_P0_mg_boundaries.push_back(values);
      const size_t index = lbs_solver.incident_P0_mg_boundaries.size()-1;

      //bid = XMIN or XMAX or YMIN ... etc
      //index is where it is on the incident_P0_mg_boundaries stack
      lbs_solver.boundary_types[bid].first =
        lbs::BoundaryType::INCIDENT_ISOTROPIC;
      lbs_solver.boundary_types[bid].second= static_cast<int>(index);

      chi::log.Log()
        << "Isotropic boundary condition for boundary " << bid
        << " loaded with " << table_len << " groups.";
    }
    else if (btype == (int)lbs::BoundaryType::REFLECTING)
    {
      lbs_solver.boundary_types[bid].first = lbs::BoundaryType::REFLECTING;
      chi::log.Log() << "Boundary " << bid << " set to Reflecting.";
    }
    else
    {
      chi::log.LogAllError()
        << "Unsupported boundary type encountered "
           "in call to " << LuaSourceInfo(L,"chiLBSSetProperty");
      chi::Exit(EXIT_FAILURE);
    }

  }
  else if (property == SCATTERING_ORDER)
  {
    LuaCheckNilValue(fname, L, 3);

    const int scattering_order = lua_tonumber(L,3);

    if (scattering_order<0)
    {
      chi::log.Log0Error()
        << "Invalid scattering order in call to "
        << "chiLBSSetProperty:SCATTERING_ORDER. "
           "Value must be > 0.";
      chi::Exit(EXIT_FAILURE);
    }

    lbs_solver.options.scattering_order = scattering_order;
  }
  else if (property == SWEEP_EAGER_LIMIT)
  {
    if (numArgs!=3)
      LuaPostArgAmountError("chiLBSSetProperty:SWEEP_EAGER_LIMIT",
                            3,numArgs);

    LuaCheckNilValue(fname, L, 3);

    const int limit = lua_tonumber(L,3);
    lbs_solver.options.sweep_eager_limit = limit;
  }
  else if (property == READ_RESTART_DATA)
  {
    if (numArgs >= 3)
    {
      LuaCheckNilValue(fname, L, 3);

      const std::string folder = lua_tostring(L,3);
      lbs_solver.options.read_restart_folder_name = std::string(folder);
      chi::log.Log() << "Restart input folder set to " << folder;
    }
    if (numArgs >= 4)
    {
      LuaCheckNilValue(fname, L, 4);

      const std::string filebase = lua_tostring(L,4);
      lbs_solver.options.read_restart_file_base = std::string(filebase);
      chi::log.Log() << "Restart input filebase set to " << filebase;
    }
    lbs_solver.options.read_restart_data = true;
  }
  else if (property == WRITE_RESTART_DATA)
  {
    if (numArgs >= 3)
    {
      LuaCheckNilValue(fname, L, 3);

      const std::string folder = lua_tostring(L,3);
      lbs_solver.options.write_restart_folder_name = std::string(folder);
      chi::log.Log() << "Restart output folder set to " << folder;
    }
    if (numArgs >= 4)
    {
      LuaCheckNilValue(fname, L, 4);

      const std::string filebase = lua_tostring(L,4);
      lbs_solver.options.write_restart_file_base = std::string(filebase);
      chi::log.Log() << "Restart output filebase set to " << filebase;
    }
    if (numArgs == 5)
    {
      LuaCheckNilValue(fname, L, 5);

      const double interval = lua_tonumber(L,5);
      lbs_solver.options.write_restart_interval = interval;
    }
    lbs_solver.options.write_restart_data = true;
  }
  else if (property == SAVE_ANGULAR_FLUX)
  {
    LuaCheckNilValue(fname, L, 3);

    const bool save_flag = lua_toboolean(L, 3);

    lbs_solver.options.save_angular_flux = save_flag;

    chi::log.Log() << "LBS option to save angular flux set to " << save_flag;
  }
  else if (property == USE_SOURCE_MOMENTS)
  {
    LuaCheckNilValue(fname, L, 3);

    const bool use_flag = lua_toboolean(L, 3);

    lbs_solver.options.use_src_moments = use_flag;

    chi::log.Log() << "LBS option to use source moments set to " << use_flag;
  }
  else if (property == VERBOSE_INNER_ITERATIONS)
  {
    LuaCheckNilValue(fname, L, 3);

    const bool flag = lua_toboolean(L, 3);

    lbs_solver.options.verbose_inner_iterations = flag;

    chi::log.Log() << "LBS option: verbose_inner_iterations set to " << flag;
  }
  else if (property == VERBOSE_OUTER_ITERATIONS)
  {
    LuaCheckNilValue(fname, L, 3);

    const bool flag = lua_toboolean(L, 3);

    lbs_solver.options.verbose_outer_iterations = flag;

    chi::log.Log() << "LBS option: verbose_outer_iterations set to " << flag;
  }
  else if (property == USE_PRECURSORS)
  {
    LuaCheckNilValue(fname, L, 3);

    const bool flag = lua_toboolean(L, 3);

    lbs_solver.options.use_precursors = flag;

    chi::log.Log() << "LBS option: use_precursors set to " << flag;
  }
  else
    throw std::logic_error(fname+": Invalid property in chiLBSSetProperty.\n");

  return 0;
}
