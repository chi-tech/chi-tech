#include "ChiLua/chi_lua.h"

#include "ChiMesh/FieldFunctionInterpolation/chi_ffinterpolation.h"

//###################################################################
/** Initialize interpolator.
 *
\param FFIHandle int Handle to the field function interpolation.

\ingroup LuaFFInterpol
\author Jan*/
int chiFFInterpolationInitialize(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  //================================================== Get handle to field function
  const size_t ffihandle = lua_tonumber(L,1);

  auto p_ffi = chi::GetStackItemPtr(chi::field_func_interpolation_stack,
                                    ffihandle, fname);

  p_ffi->Initialize();
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
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  //================================================== Get handle to field function
  const size_t ffihandle = lua_tonumber(L,1);

  auto p_ffi = chi::GetStackItemPtr(chi::field_func_interpolation_stack,
                                    ffihandle, fname);

  p_ffi->Execute();
  return 0;
}

