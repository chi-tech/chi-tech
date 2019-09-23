#include "ChiLua/chi_lua.h"

#include "../lbs_linear_boltzman_solver.h"
#include "ChiPhysics/chi_physics.h"
#include "ChiMath/chi_math.h"

#include "ChiMath/SpatialDiscretization/PiecewiseLinear/pwl.h"

extern ChiPhysics chi_physics_handler;
extern ChiMath    chi_math_handler;

#define DISCRETIZATION_METHOD 1
  #define PWLD1D 4
  #define PWLD2D 5
  #define PWLD3D 6

#define PARTITION_METHOD 2
  #define SERIAL 1
  #define FROM_SURFACE 2

#define BOUNDARY_CONDITION 3
  #define XMAX 31
  #define XMIN 32
  #define YMAX 33
  #define YMIN 34
  #define ZMAX 35
  #define ZMIN 36
//Already defined in chi_nptransport.h but just here for convenience
//    #define VACUUM             301
//    #define INCIDENT_ISOTROPIC 302
#define GROUPSET_TOLERANCE     102
#define GROUPSET_MAXITERATIONS 103
#define GROUPSET_GMRESRESTART_INTVL 104
#define GROUPSET_SUBSETS 105
#define GROUPSET_WGDSA   106
#define GROUPSET_TGDSA   107
#define GROUPSET_WGDSA_MAXITERATIONS 108
#define GROUPSET_TGDSA_MAXITERATIONS 109
#define GROUPSET_WGDSA_TOLERANCE 110
#define GROUPSET_TGDSA_TOLERANCE 111

#define SCATTERING_ORDER 4

#define SWEEP_EAGER_LIMIT 5


//###################################################################
/**Set LBS property.
\param SolverIndex int Handle to the solver for which the set is to be created.
\param PropertyIndex int Code for a specific property.


##_

###PropertyIndex\n
DISCRETIZATION_METHOD\n
 Discretization method.\n\n

PARTITION_METHOD\n
 Multi-processor partitioning method.\n\n

BOUNDARY_CONDITION\n
 Boundary condition type. See BoundaryIdentify.\n\n

GROUPSET_ITERATIVEMETHOD\n
 Iterative method to be used by a given groupset. Expects to be followed
 by the groupset handle and IterativeMethod.\n\n

GROUPSET_TOLERANCE\n
 Either Residual or pointwise tolerance to use for iterative method. Expects
 to be followed by the groupset handle and a float. Default is 1.0e-6\n\n

GROUPSET_MAXITERATIONS\n
 Maximum iterations for groupset solve. Exepects to be followed by the groupset
 handle and an integer. Default is 1000.\n\n

GROUPSET_GMRESRESTART_INTVL\n
 Number of GMRES iterations to accumulate Krylov vectors before restarting
 the accumulation. Default 10.\n\n

GROUPSET_SUBSETS\n
 Number of subsets to use for groupset. Default 1. Expects to be followed by
 groupset handle and integer amount of subsets.\n\n

GROUPSET_WGDSA\n
 Expects to be followed by a boolean flag. If true then Within-Group
 Diffusion Synthetic Acceleration will be applied. Default false.\n\n

GROUPSET_TGDSA\n
 Expects to be followed by a boolean flag. If true then Two-Grid
 Diffusion Synthetic Acceleration will be applied. Default false.\n\n

SCATTERING_ORDER\n
 Defines the level of harmonic expansion for the scattering source.Default 1.
 Expects to be followed by an integer.\n\n

SWEEP_EAGER_LIMIT\n
 The eager limit to be used in message size during sweep initialization.
 This expects to be followed by a size in bytes (Max 64,0000).Default 32,000.
 See note below.\n\n


###Discretization methods
 PWLD2D = Piecewise Linear Finite Element 2D.\n
 PWLD3D = Piecewise Linear Finite Element 3D.

###Partitioning methods
 SERIAL = No multi-processing.\n
 FROM_SURFACE = Same partitioning as used on Surface mesh.

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
VACUUM\n
Specifies a vaccuum boundary condition. It is not followed by any value.\n
\n
\n
INCIDENT_ISOTROPIC\n
Incident isotropic flux. This argument needs to be followed by a lua table
index 1 to G where G is the amount of energy groups. Note internally this
is mapped as 0 to G-1.

###IterativeMethod
NPT_CLASSICRICHARDSON\n
Classic richardson iteration or otherwise known as source iteration.\n
\n
NPT_GMRES\n
Generalized Minimized Residual. Very robust method for converging solutions.\n

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

\ingroup LuaNPT*/
int chiLBSSetProperty(lua_State *L)
{
  int numArgs = lua_gettop(L);
  int solver_index = lua_tonumber(L,1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      fprintf(stderr,"ERROR: Incorrect solver-type"
                     "in chiLBSSetProperty\n");
      exit(EXIT_FAILURE);
    }
  }
  catch(const std::out_of_range& o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiLBSSetProperty\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Get property index
  int property = lua_tonumber(L,2);

  //============================================= Handle properties
  if (property == DISCRETIZATION_METHOD)
  {
    int method = lua_tonumber(L,3);
    if (method == PWLD1D)
    {
      SpatialDiscretization_PWL* discretization = new SpatialDiscretization_PWL(1);
      solver->discretization = discretization;
    }
    else if (method == PWLD2D)
    {
      SpatialDiscretization_PWL* discretization = new SpatialDiscretization_PWL(2);
      solver->discretization = discretization;
    }
    else if (method == PWLD3D)
    {
      SpatialDiscretization_PWL* discretization = new SpatialDiscretization_PWL(3);
      solver->discretization = discretization;
    }
    else
    {
      std::cerr << "Invalid option for Discretization method in "
                   "chiLBSSetProperty.\n";
      exit(EXIT_FAILURE);
    }
  }
  else if (property == PARTITION_METHOD)
  {
    int method = lua_tonumber(L,3);
    solver->options.partition_method = method;
    //printf("Partition method set to %d\n",method);
  }
  else if (property == BOUNDARY_CONDITION)
  {
    if (numArgs<4)
      LuaPostArgAmountError("chiLBSSetProperty",4,numArgs);

    int bident = lua_tonumber(L,3);
    int btype  = lua_tonumber(L,4);

    if (!((bident>=XMAX) && (bident<=ZMIN)))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Unknown boundary identifier encountered "
           "in call to chiLBSSetProperty";
      exit(EXIT_FAILURE);
    }

    int bid = bident - 31;

    if (btype == VACUUM)
    {
      solver->boundary_types[bid].first = VACUUM;
    }
    else if (btype == INCIDENT_ISOTROPIC)
    {
      if (numArgs!=5)
        LuaPostArgAmountError("chiLBSSetProperty",5,numArgs);

      if (solver->groups.size() == 0)
      {
        chi_log.Log(LOG_0ERROR)
          << "In call to chiLBSSetProperty, setting "
          << "incident isotropic flux boundary type: Number of solver groups"
          << " is zero. Boundary fluxes can only be set after group structure"
          << " has been defined.";
        exit(EXIT_FAILURE);
      }

      if (!lua_istable(L,5))
      {
        chi_log.Log(LOG_ALLERROR)
          << "In call to chiLBSSetProperty, setting "
          << "incident isotropic flux boundary type,"
          << " argument 5 should be a lua table and was detected as"
             " not being one.";
        exit(EXIT_FAILURE);
      }

      int table_len = lua_rawlen(L,5);
      std::vector<double> values(table_len,0.0);
      for (int g=0; g<table_len; g++)
      {
        lua_pushnumber(L,g+1);
        lua_gettable(L,5);
        values[g] = lua_tonumber(L,-1);
        lua_pop(L,1);
//        chi_log.Log(LOG_0) << g << " " << values[g];
      }

      if (table_len != solver->groups.size())
      {
        chi_log.Log(LOG_0ERROR)
          << "In call to chiLBSSetProperty, setting "
          << "incident isotropic flux boundary type: "
          << "Number of groups in boundary flux specification is "
          << table_len << " but solver has a total of "
          << solver->groups.size() << " groups. These two must be equal.";
        exit(EXIT_FAILURE);
      }

      solver->incident_P0_mg_boundaries.push_back(values);
      int index = solver->incident_P0_mg_boundaries.size()-1;

      //bid = XMIN or XMAX or YMIN ... etc
      //index is where it is on the incident_P0_mg_boundaries stack
      solver->boundary_types[bid].first = INCIDENT_ISOTROPIC;
      solver->boundary_types[bid].second= index;

      chi_log.Log(LOG_0)
        << "Isotropic boundary condition for boundary " << bid
        << " loaded with " << table_len << " groups.";
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Unsupported boundary type encountered "
           "in call to chiLBSSetProperty";
      exit(EXIT_FAILURE);
    }

  }
  else if (property == GROUPSET_ITERATIVEMETHOD)
  {
    if (numArgs!=4)
      LuaPostArgAmountError("chiLBSSetProperty:GROUPSET_ITERATIVEMETHOD",
                            4,numArgs);

    int groupset_num = lua_tonumber(L,3);
    LBSGroupset* groupset;
    try {
      groupset = solver->group_sets.at(groupset_num);
    }
    catch(const std::out_of_range& o){
      chi_log.Log(LOG_0ERROR)
        << "Invalid Groupset handle supplied in call to "
        << "chiLBSSetProperty:GROUPSET_ITERATIVEMETHOD.";
      exit(EXIT_FAILURE);
    }

    int iter_method = lua_tonumber(L,4);

    if (iter_method == NPT_CLASSICRICHARDSON)
    {
      groupset->iterative_method = NPT_CLASSICRICHARDSON;
    }
    else if (iter_method == NPT_GMRES)
    {
      groupset->iterative_method = NPT_GMRES;
    }
    else
    {
      chi_log.Log(LOG_0ERROR)
        << "Unsupported iterative method specified in call to "
        << "chiLBSSetProperty:GROUPSET_ITERATIVEMETHOD.";
      exit(EXIT_FAILURE);
    }

  }
  else if (property == GROUPSET_TOLERANCE)
  {
    if (numArgs!=4)
      LuaPostArgAmountError("chiLBSSetProperty:GROUPSET_TOLERANCE",
                            4,numArgs);

    int groupset_num = lua_tonumber(L,3);
    LBSGroupset* groupset;
    try {
      groupset = solver->group_sets.at(groupset_num);
    }
    catch(const std::out_of_range& o){
      chi_log.Log(LOG_0ERROR)
        << "Invalid Groupset handle supplied in call to "
        << "chiLBSSetProperty:GROUPSET_TOLERANCE.";
      exit(EXIT_FAILURE);
    }

    groupset->residual_tolerance = lua_tonumber(L,4);
  }
  else if (property == GROUPSET_MAXITERATIONS)
  {
    if (numArgs!=4)
      LuaPostArgAmountError("chiLBSSetProperty:GROUPSET_MAXITERATIONS",
                            4,numArgs);

    int groupset_num = lua_tonumber(L,3);
    LBSGroupset* groupset;
    try {
      groupset = solver->group_sets.at(groupset_num);
    }
    catch(const std::out_of_range& o){
      chi_log.Log(LOG_0ERROR)
        << "Invalid Groupset handle supplied in call to "
        << "chiLBSSetProperty:GROUPSET_MAXITERATIONS.";
      exit(EXIT_FAILURE);
    }

    groupset->max_iterations = lua_tonumber(L,4);
  }
  else if (property == GROUPSET_GMRESRESTART_INTVL)
  {
    if (numArgs!=4)
      LuaPostArgAmountError("chiLBSSetProperty:GROUPSET_GMRESRESTART_INTVL",
                            4,numArgs);

    int groupset_num = lua_tonumber(L,3);
    LBSGroupset* groupset;
    try {
      groupset = solver->group_sets.at(groupset_num);
    }
    catch(const std::out_of_range& o){
      chi_log.Log(LOG_0ERROR)
        << "Invalid Groupset handle supplied in call to "
        << "chiLBSSetProperty:GROUPSET_GMRESRESTART_INTVL.";
      exit(EXIT_FAILURE);
    }

    groupset->gmres_restart_intvl = lua_tonumber(L,4);
  }
  else if (property == GROUPSET_SUBSETS)
  {
    if (numArgs!=4)
      LuaPostArgAmountError("chiLBSSetProperty:GROUPSET_SUBSETS",
                            4,numArgs);

    int groupset_num = lua_tonumber(L,3);
    LBSGroupset* groupset;
    try {
      groupset = solver->group_sets.at(groupset_num);
    }
    catch(const std::out_of_range& o){
      chi_log.Log(LOG_0ERROR)
        << "Invalid Groupset handle supplied in call to "
        << "chiLBSSetProperty:GROUPSET_SUBSETS.";
      exit(EXIT_FAILURE);
    }

    groupset->master_num_grp_subsets = lua_tonumber(L,4);
  }
  else if (property == GROUPSET_WGDSA)
  {
    if (numArgs!=4)
      LuaPostArgAmountError("chiLBSSetProperty:GROUPSET_WGDSA",
                            4,numArgs);

    int groupset_num = lua_tonumber(L,3);
    LBSGroupset* groupset;
    try {
      groupset = solver->group_sets.at(groupset_num);
    }
    catch(const std::out_of_range& o){
      chi_log.Log(LOG_0ERROR)
        << "Invalid Groupset handle supplied in call to "
        << "chiLBSSetProperty:GROUPSET_WGDSA.";
      exit(EXIT_FAILURE);
    }

    if (!lua_isboolean(L,4))
    {
      chi_log.Log(LOG_0ERROR)
        << "Invalid datatype supplied in call to "
        << "chiLBSSetProperty:GROUPSET_WGDSA. Boolean expected.";
      exit(EXIT_FAILURE);
    }

    groupset->apply_wgdsa = lua_toboolean(L,4);
  }
  else if (property == GROUPSET_TGDSA)
  {
    if (numArgs!=4)
      LuaPostArgAmountError("chiLBSSetProperty:GROUPSET_TGDSA",
                            4,numArgs);

    int groupset_num = lua_tonumber(L,3);
    LBSGroupset* groupset;
    try {
      groupset = solver->group_sets.at(groupset_num);
    }
    catch(const std::out_of_range& o){
      chi_log.Log(LOG_0ERROR)
        << "Invalid Groupset handle supplied in call to "
        << "chiLBSSetProperty:GROUPSET_TGDSA.";
      exit(EXIT_FAILURE);
    }

    if (!lua_isboolean(L,4))
    {
      chi_log.Log(LOG_0ERROR)
        << "Invalid datatype supplied in call to "
        << "chiLBSSetProperty:GROUPSET_TGDSA. Boolean expected.";
      exit(EXIT_FAILURE);
    }

    groupset->apply_tgdsa = lua_toboolean(L,4);
  }
  else if (property == GROUPSET_WGDSA_MAXITERATIONS)
  {
    if (numArgs!=4)
      LuaPostArgAmountError("chiLBSSetProperty:GROUPSET_WGDSA_MAXITERATIONS",
                            4,numArgs);

    int groupset_num = lua_tonumber(L,3);
    LBSGroupset* groupset;
    try {
      groupset = solver->group_sets.at(groupset_num);
    }
    catch(const std::out_of_range& o){
      chi_log.Log(LOG_0ERROR)
        << "Invalid Groupset handle supplied in call to "
        << "chiLBSSetProperty:GROUPSET_WGDSA_MAXITERATIONS.";
      exit(EXIT_FAILURE);
    }

    groupset->wgdsa_max_iters = lua_tonumber(L,4);
  }
  else if (property == GROUPSET_TGDSA_MAXITERATIONS)
  {
    if (numArgs!=4)
      LuaPostArgAmountError("chiLBSSetProperty:GROUPSET_TGDSA_MAXITERATIONS",
                            4,numArgs);

    int groupset_num = lua_tonumber(L,3);
    LBSGroupset* groupset;
    try {
      groupset = solver->group_sets.at(groupset_num);
    }
    catch(const std::out_of_range& o){
      chi_log.Log(LOG_0ERROR)
        << "Invalid Groupset handle supplied in call to "
        << "chiLBSSetProperty:GROUPSET_TGDSA_MAXITERATIONS.";
      exit(EXIT_FAILURE);
    }

    groupset->tgdsa_max_iters = lua_tonumber(L,4);
  }
  else if (property == GROUPSET_WGDSA_TOLERANCE)
  {
    if (numArgs!=4)
      LuaPostArgAmountError("chiLBSSetProperty:GROUPSET_WGDSA_TOLERANCE",
                            4,numArgs);

    int groupset_num = lua_tonumber(L,3);
    LBSGroupset* groupset;
    try {
      groupset = solver->group_sets.at(groupset_num);
    }
    catch(const std::out_of_range& o){
      chi_log.Log(LOG_0ERROR)
        << "Invalid Groupset handle supplied in call to "
        << "chiLBSSetProperty:GROUPSET_WGDSA_TOLERANCE.";
      exit(EXIT_FAILURE);
    }

    groupset->wgdsa_tol = lua_tonumber(L,4);
  }
  else if (property == GROUPSET_TGDSA_TOLERANCE)
  {
    if (numArgs!=4)
      LuaPostArgAmountError("chiLBSSetProperty:GROUPSET_TGDSA_TOLERANCE",
                            4,numArgs);

    int groupset_num = lua_tonumber(L,3);
    LBSGroupset* groupset;
    try {
      groupset = solver->group_sets.at(groupset_num);
    }
    catch(const std::out_of_range& o){
      chi_log.Log(LOG_0ERROR)
        << "Invalid Groupset handle supplied in call to "
        << "chiLBSSetProperty:GROUPSET_TGDSA_TOLERANCE.";
      exit(EXIT_FAILURE);
    }

    groupset->tgdsa_tol = lua_tonumber(L,4);
  }
  else if (property == SCATTERING_ORDER)
  {
    int scat_order = lua_tonumber(L,3);

    if (scat_order<0)
    {
      chi_log.Log(LOG_0ERROR)
        << "Invalid scattering order in call to "
        << "chiLBSSetProperty:SCATTERING_ORDER. "
           "Value must be > 0.";
      exit(EXIT_FAILURE);
    }

    solver->options.scattering_order = scat_order;
  }
  else if (property == SWEEP_EAGER_LIMIT)
  {
    if (numArgs!=3)
      LuaPostArgAmountError("chiLBSSetProperty:SWEEP_EAGER_LIMIT",
                            3,numArgs);

    int limit = lua_tonumber(L,3);
    if (limit<=64000)
    {
      solver->options.sweep_eager_limit = limit;
    }
  }
  else
  {
    std::cerr << "Invalid property in chiLBSSetProperty.\n";
    exit(EXIT_FAILURE);
  }

  return 0;
}
