#include "ChiLua/chi_lua.h"
#include "ChiMesh/FieldFunctionInterpolation/Slice/chi_ffinter_slice.h"
#include "ChiMesh/FieldFunctionInterpolation/Line/chi_ffinter_line.h"
#include "ChiMesh/FieldFunctionInterpolation/Volume/chi_ffinter_volume.h"

#include "chi_runtime.h"

#include "chi_runtime.h"
#include "chi_log.h"
;


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
  const std::string fname = __FUNCTION__;

  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiFFInterpolationGetValue",1,num_args);

  double value = 0.0;

  //================================================== Get handle to field function
  const size_t ffihandle = lua_tonumber(L,1);

  auto p_ffi = chi::GetStackItemPtr(chi::field_func_interpolation_stack,
                                    ffihandle, fname);

  if (typeid(*p_ffi) == typeid(chi_mesh::FieldFunctionInterpolationVolume))
  {
    auto& cur_ffi_volume = (chi_mesh::FieldFunctionInterpolationVolume&)*p_ffi;
    value = cur_ffi_volume.op_value;

    lua_pushnumber(L,value);
    return 1;
  }
  else if (typeid(*p_ffi) == typeid(chi_mesh::FieldFunctionInterpolationLine))
  {
    auto& cur_ffi_line = (chi_mesh::FieldFunctionInterpolationLine&)*p_ffi;

    lua_newtable(L);

    for (int ff=0; ff<cur_ffi_line.field_functions.size(); ff++)
    {
      lua_pushnumber(L,ff+1);

      lua_newtable(L);
      auto ff_ctx = cur_ffi_line.ff_contexts[ff];

      for (int p=0; p<cur_ffi_line.interpolation_points.size(); p++)
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
    chi::log.Log0Warning()
      << "chiFFInterpolationGetValue is currently only supported for "
      << " VOLUME interpolator types.";
  }

  return 0;
}
