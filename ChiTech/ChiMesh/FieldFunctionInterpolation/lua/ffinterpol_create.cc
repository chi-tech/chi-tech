#include "ChiLua/chi_lua.h"
#include "ChiMesh/FieldFunctionInterpolation/Slice/chi_ffinter_slice.h"
#include "ChiMesh/FieldFunctionInterpolation/Line/chi_ffinter_line.h"
#include "ChiMesh/FieldFunctionInterpolation/Volume/chi_ffinter_volume.h"

#include "chi_runtime.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;



//#############################################################################
/** Creates a new field function interpolation.
 *
\param FFITypeIndex int Type of field function interpolation.

##_

###FFITypeIndex\n
SLICE           = Two dimensional slice of the mesh. \n
LINE            = Line defined by two points.\n
VOLUME          = Volume either referring to the entire volume or that of a
                  logical volume assigned to the interpolator.\n

\return Handle int Handle to the created interpolation.
\ingroup LuaFFInterpol
\author Jan*/
int chiFFInterpolationCreate(lua_State *L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  //================================================== Process types
  int ffitype = lua_tonumber(L,1);
  if (ffitype == FFI_SLICE)                         //SLICE
  {
    auto new_ffi = new chi_mesh::FieldFunctionInterpolationSlice;

    chi::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = chi::field_func_interpolation_stack.size()-1;
    chi_log.Log(LOG_ALLVERBOSE_2)
    << "Created slice Field Function Interpolation";
    lua_pushnumber(L,static_cast<lua_Number>(index));
    return 1;
  }
  else if (ffitype == FFI_LINE)
  {
    auto new_ffi = new chi_mesh::FieldFunctionInterpolationLine;

    chi::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = chi::field_func_interpolation_stack.size()-1;
    chi_log.Log(LOG_ALLVERBOSE_2)
      << "Created line Field Function Interpolation";
    lua_pushnumber(L,static_cast<lua_Number>(index));
    return 1;
  }
  else if (ffitype == FFI_VOLUME)
  {
    auto new_ffi = new chi_mesh::FieldFunctionInterpolationVolume;

    chi::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = chi::field_func_interpolation_stack.size()-1;
    chi_log.Log(LOG_ALLVERBOSE_2)
      << "Created Volume Field Function Interpolation";
    lua_pushnumber(L,static_cast<lua_Number>(index));
    return 1;
  }
  else                                              //Fall back
  {
    chi_log.Log(LOG_ALLERROR)
    << "Invalid FFITypeIndex used in chiFFInterpolationCreate.";
    exit(EXIT_FAILURE);
  }
}