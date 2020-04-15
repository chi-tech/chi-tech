#include "../../../ChiLua/chi_lua.h"
#include "../../MeshHandler/chi_meshhandler.h"
#include "../../FieldFunctionInterpolation/Slice/chi_ffinter_slice.h"
#include "../../FieldFunctionInterpolation/Line/chi_ffinter_line.h"
#include "../../FieldFunctionInterpolation/Volume/chi_ffinter_volume.h"
#include <ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h>
#include "../../../ChiPhysics/chi_physics.h"

#include <chi_log.h>

extern ChiLog chi_log;
extern ChiPhysics chi_physics_handler;

//#############################################################################
/** Gets the value(s) associated with an interpolation provided the
 * interpolation type has an associated value.
 *
\param FFIHandle int Handle to the field function interpolation.

###Note:
Currently only the Volume interpolation supports obtaining a value.

\ingroup LuaFFInterpol
\author Jan*/
int chiFFInterpolationGetValue(lua_State *L)
{
  double value = 0.0;

  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiFFInterpolationGetValue",1,num_args);

  //================================================== Get handle to field function
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();
  int ffihandle = lua_tonumber(L,1);
  chi_mesh::FieldFunctionInterpolation* cur_ffi;
  try {
    cur_ffi = cur_hndlr->ffinterpolation_stack.at(ffihandle);
  }
  catch(const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid ffi handle in chiFFInterpolationGetValue.";
    exit(EXIT_FAILURE);
  }

  if (typeid(*cur_ffi) == typeid(chi_mesh::FieldFunctionInterpolationVolume))
  {
    auto cur_ffi_volume = (chi_mesh::FieldFunctionInterpolationVolume*)cur_ffi;
    value = cur_ffi_volume->op_value;

    lua_pushnumber(L,value);
    return 1;
  }
  else if (typeid(*cur_ffi) == typeid(chi_mesh::FieldFunctionInterpolationLine))
  {
    auto cur_ffi_line = (chi_mesh::FieldFunctionInterpolationLine*)cur_ffi;

    lua_newtable(L);

    for (int ff=0; ff<cur_ffi_line->field_functions.size(); ff++)
    {
      lua_pushnumber(L,ff+1);

      lua_newtable(L);
      auto ff_ctx = cur_ffi_line->ff_contexts[ff];

      for (int p=0; p<cur_ffi_line->interpolation_points.size(); p++)
      {
        lua_pushnumber(L,p+1);
        lua_pushnumber(L,ff_ctx->interpolation_points_values[p]);
        lua_settable(L,-3);
      }

      lua_settable(L,-3);
    }

    return 1;
  }
  else
  {
    chi_log.Log(LOG_0WARNING)
      << "chiFFInterpolationGetValue is currently only supported for "
      << " VOLUME interpolator types.";
  }

  return 0;
}
