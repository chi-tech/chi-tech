#include "../../../CHI_LUA/chi_lua.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../../CHI_PHYSICS/chi_physics.h"
#include "../chi_ffinterpolation.h"

#include <chi_log.h>

extern CHI_LOG chi_log;
extern CHI_PHYSICS chi_physics_handler;

//###################################################################
/** Initialize interpolator.
 *
\param FFIHandle int Handle to the field function interpolation.

\ingroup LuaFFInterpol
\author Jan*/
int chiFFInterpolationInitialize(lua_State* L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //================================================== Get handle to field function
  int ffihandle = lua_tonumber(L,1);
  chi_mesh::FieldFunctionInterpolation* cur_ffi;
  try {
    cur_ffi = cur_hndlr->ffinterpolation_stack.at(ffihandle);
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid ffi handle in chiFFInterpolationSetProperty.";
    exit(EXIT_FAILURE);
  }

  cur_ffi->Initialize();
  return 0;
}

//###################################################################
/** Execute interpolator.
 *
\param FFIHandle int Handle to the field function interpolation.

\ingroup LuaFFInterpol
\author Jan*/
int chiFFInterpolationExecute(lua_State* L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //================================================== Get handle to field function
  int ffihandle = lua_tonumber(L,1);
  chi_mesh::FieldFunctionInterpolation* cur_ffi;
  try {
    cur_ffi = cur_hndlr->ffinterpolation_stack.at(ffihandle);
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid ffi handle in chiFFInterpolationSetProperty.";
    exit(EXIT_FAILURE);
  }

  cur_ffi->Execute();
  return 0;
}

