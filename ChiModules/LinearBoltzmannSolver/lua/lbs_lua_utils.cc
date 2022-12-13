#include "lbs_lua_utils.h"

#include "chi_runtime.h"

lbs::SteadySolver& lbs::lua_utils::
  GetSolverByHandle(int handle, const std::string& calling_function_name)
{
  std::shared_ptr<lbs::SteadySolver> lbs_solver;
  try{

    lbs_solver = std::dynamic_pointer_cast<lbs::SteadySolver>(
      chi::solver_stack.at(handle));

    if (not lbs_solver)
      throw std::logic_error(calling_function_name +
      ": Invalid solver at given handle (" +
      std::to_string(handle) + "). "
      "The solver is not of type LinearBoltzmann::Solver.");
  }//try
  catch(const std::out_of_range& o) {
    throw std::logic_error(calling_function_name + ": Invalid solver-handle (" +
                           std::to_string(handle) + ").");
  }

  return *lbs_solver;
}

#define LUA_FMACRO1(x) lua_register(L, #x, x)
#define LUA_CMACRO1(x,y) \
        lua_pushnumber(L, y); \
        lua_setglobal(L, #x)

#define LUA_CTABLE1(x) \
        lua_newtable(L); \
        lua_setglobal(L, #x)

#define LUA_CADDCONST_VALUE_TO_TABLE1(const_name,const_value,namespace_name) \
        lua_getglobal(L,#namespace_name); \
        lua_pushstring(L,#const_name); \
        lua_pushnumber(L,const_value); \
        lua_settable(L,-3); \
        lua_pop(L,1)

void lbs::lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiLBSCreateSolver);

  LUA_FMACRO1(chiLBSSetProperty);
    LUA_CMACRO1(DISCRETIZATION_METHOD,   1);
    LUA_CMACRO1(PWLD, 3);
    LUA_CMACRO1(PWLD1D, 4);
    LUA_CMACRO1(PWLD2D, 5);
    LUA_CMACRO1(PWLD3D, 6);
    LUA_CMACRO1(PARTITION_METHOD,   2);
    LUA_CMACRO1(SERIAL,   1);
    LUA_CMACRO1(FROM_SURFACE,   2);
    LUA_CMACRO1(BOUNDARY_CONDITION,   3);
    LUA_CMACRO1(XMAX,   31);
    LUA_CMACRO1(XMIN,   32);
    LUA_CMACRO1(YMAX,   33);
    LUA_CMACRO1(YMIN,   34);
    LUA_CMACRO1(ZMAX,   35);
    LUA_CMACRO1(ZMIN,   36);
    LUA_CMACRO1(SCATTERING_ORDER,   4);
    LUA_CMACRO1(SWEEP_EAGER_LIMIT,   5);
    LUA_CMACRO1(READ_RESTART_DATA,   6);
    LUA_CMACRO1(WRITE_RESTART_DATA,  7);
    LUA_CMACRO1(SAVE_ANGULAR_FLUX, 8);
    LUA_CMACRO1(USE_SOURCE_MOMENTS, 9);
    LUA_CMACRO1(VERBOSE_INNER_ITERATIONS, 10);
    LUA_CMACRO1(VERBOSE_OUTER_ITERATIONS, 11);
    LUA_CMACRO1(USE_PRECURSORS, 12);


  LUA_CTABLE1(LBSProperty);
  LUA_CADDCONST_VALUE_TO_TABLE1(DISCRETIZATION_METHOD, 1, LBSProperty);
  LUA_CADDCONST_VALUE_TO_TABLE1(BOUNDARY_CONDITION,    3, LBSProperty);
  LUA_CADDCONST_VALUE_TO_TABLE1(SCATTERING_ORDER,      4, LBSProperty);
  LUA_CADDCONST_VALUE_TO_TABLE1(SWEEP_EAGER_LIMIT,     5, LBSProperty);
  LUA_CADDCONST_VALUE_TO_TABLE1(READ_RESTART_DATA,     6, LBSProperty);
  LUA_CADDCONST_VALUE_TO_TABLE1(WRITE_RESTART_DATA,    7, LBSProperty);
  LUA_CADDCONST_VALUE_TO_TABLE1(SAVE_ANGULAR_FLUX,     8, LBSProperty);

  LUA_CTABLE1(LBSSpatialDiscretizations);
  LUA_CADDCONST_VALUE_TO_TABLE1(PWLD, 3, LBSSpatialDiscretizations);

  LUA_CTABLE1(LBSBoundaryID);
  LUA_CADDCONST_VALUE_TO_TABLE1(XMAX,   31, LBSBoundaryID);
  LUA_CADDCONST_VALUE_TO_TABLE1(XMIN,   32, LBSBoundaryID);
  LUA_CADDCONST_VALUE_TO_TABLE1(YMAX,   33, LBSBoundaryID);
  LUA_CADDCONST_VALUE_TO_TABLE1(YMIN,   34, LBSBoundaryID);
  LUA_CADDCONST_VALUE_TO_TABLE1(ZMAX,   35, LBSBoundaryID);
  LUA_CADDCONST_VALUE_TO_TABLE1(ZMIN,   36, LBSBoundaryID);

  LUA_CTABLE1(LBSBoundaryTypes);
  LUA_CADDCONST_VALUE_TO_TABLE1(VACUUM            ,1,LBSBoundaryTypes);
  LUA_CADDCONST_VALUE_TO_TABLE1(INCIDENT_ISOTROPIC,2,LBSBoundaryTypes);
  LUA_CADDCONST_VALUE_TO_TABLE1(REFLECTING        ,3,LBSBoundaryTypes);

  LUA_CMACRO1(GROUPSET_ITERATIVEMETHOD    , 101);
  LUA_CMACRO1(NPT_CLASSICRICHARDSON       , 1);
  LUA_CMACRO1(NPT_CLASSICRICHARDSON_CYCLES, 2);
  LUA_CMACRO1(NPT_GMRES                   , 3);
  LUA_CMACRO1(NPT_GMRES_CYCLES            , 4);

  LUA_CMACRO1(KRYLOV_RICHARDSON           , 5);
  LUA_CMACRO1(KRYLOV_RICHARDSON_CYCLES    , 6);
  LUA_CMACRO1(KRYLOV_GMRES                , 7);
  LUA_CMACRO1(KRYLOV_GMRES_CYCLES         , 8);
  LUA_CMACRO1(KRYLOV_BICGSTAB             , 9);
  LUA_CMACRO1(KRYLOV_BICGSTAB_CYCLES      , 10);

  LUA_CMACRO1(GROUPSET_TOLERANCE          , 102);
  LUA_CMACRO1(GROUPSET_MAXITERATIONS      , 103);
  LUA_CMACRO1(GROUPSET_GMRESRESTART_INTVL , 104);
  LUA_CMACRO1(GROUPSET_SUBSETS            , 105);
  LUA_CMACRO1(GROUPSET_WGDSA              , 106);
  LUA_CMACRO1(GROUPSET_TGDSA              , 107);
  LUA_CMACRO1(GROUPSET_WGDSA_MAXITERATIONS, 108);
  LUA_CMACRO1(GROUPSET_TGDSA_MAXITERATIONS, 109);
  LUA_CMACRO1(GROUPSET_WGDSA_TOLERANCE    , 110);
  LUA_CMACRO1(GROUPSET_TGDSA_TOLERANCE    , 111);

  LUA_FMACRO1(chiLBSInitialize);
  LUA_FMACRO1(chiLBSExecute);
  LUA_FMACRO1(chiLBSGetFieldFunctionList);
  LUA_FMACRO1(chiLBSGetScalarFieldFunctionList);
  LUA_FMACRO1(chiLBSWriteGroupsetAngularFlux);
  LUA_FMACRO1(chiLBSReadGroupsetAngularFlux);
  LUA_FMACRO1(chiLBSWriteFluxMoments);
  LUA_FMACRO1(chiLBSCreateAndWriteSourceMoments);
  LUA_FMACRO1(chiLBSReadFluxMomentsAndMakeSourceMoments);
  LUA_FMACRO1(chiLBSReadSourceMoments);
  LUA_FMACRO1(chiLBSReadFluxMoments);
  LUA_FMACRO1(chiLBSComputeBalance);

  //=================================== Groupset manipulation
  LUA_CTABLE1(LBSGroupset);
  LUA_FMACRO1(chiLBSCreateGroup);
  LUA_FMACRO1(chiLBSCreateGroupset);
  LUA_FMACRO1(chiLBSGroupsetAddGroups);
  LUA_FMACRO1(chiLBSGroupsetSetQuadrature);
  LUA_FMACRO1(chiLBSGroupsetSetAngleAggregationType);
  LUA_CADDCONST_VALUE_TO_TABLE1(ANGLE_AGG_SINGLE,   1,LBSGroupset);
  LUA_CADDCONST_VALUE_TO_TABLE1(ANGLE_AGG_POLAR,    2,LBSGroupset);
  LUA_CADDCONST_VALUE_TO_TABLE1(ANGLE_AGG_AZIMUTHAL,3,LBSGroupset);
  LUA_FMACRO1(chiLBSGroupsetSetAngleAggDiv);
  LUA_FMACRO1(chiLBSGroupsetSetGroupSubsets);
  LUA_FMACRO1(chiLBSGroupsetSetIterativeMethod);
  LUA_FMACRO1(chiLBSGroupsetSetResidualTolerance);
  LUA_FMACRO1(chiLBSGroupsetSetMaxIterations);
  LUA_FMACRO1(chiLBSGroupsetSetGMRESRestartIntvl);
  LUA_FMACRO1(chiLBSGroupsetSetEnableSweepLog);
  LUA_FMACRO1(chiLBSGroupsetSetWGDSA);
  LUA_FMACRO1(chiLBSGroupsetSetTGDSA);

  //=================================== Point source
  LUA_FMACRO1(chiLBSAddPointSource);
  LUA_FMACRO1(chiLBSClearPointSources);
  LUA_FMACRO1(chiLBSInitializePointSources);
}
