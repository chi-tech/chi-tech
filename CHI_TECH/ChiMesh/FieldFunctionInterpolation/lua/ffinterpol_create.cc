#include "../../../ChiLua/chi_lua.h"
#include "../../MeshHandler/chi_meshhandler.h"
#include "../../FieldFunctionInterpolation/Slice/chi_ffinter_slice.h"
#include "../../FieldFunctionInterpolation/Line/chi_ffinter_line.h"
#include "../../FieldFunctionInterpolation/Volume/chi_ffinter_volume.h"
#include "../../../ChiPhysics/chi_physics.h"

/** \defgroup LuaFFInterpol Field Function Interpolation
 * \ingroup LuaMesh
*/

#include <chi_log.h>

extern ChiLog chi_log;
extern ChiPhysics chi_physics_handler;



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
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //================================================== Process types
  int ffitype = lua_tonumber(L,1);
  if (ffitype == FFI_SLICE)                         //SLICE
  {
    chi_mesh::FieldFunctionInterpolationSlice* new_ffi =
      new chi_mesh::FieldFunctionInterpolationSlice;

    cur_hndlr->ffinterpolation_stack.push_back(new_ffi);
    int index = cur_hndlr->ffinterpolation_stack.size()-1;
    chi_log.Log(LOG_ALLVERBOSE_2)
    << "Created slice Field Function Interpolation";
    lua_pushnumber(L,index);
    return 1;
  }
  else if (ffitype == FFI_LINE)
  {
    chi_mesh::FieldFunctionInterpolationLine* new_ffi =
      new chi_mesh::FieldFunctionInterpolationLine;

    cur_hndlr->ffinterpolation_stack.push_back(new_ffi);
    int index = cur_hndlr->ffinterpolation_stack.size()-1;
    chi_log.Log(LOG_ALLVERBOSE_2)
      << "Created line Field Function Interpolation";
    lua_pushnumber(L,index);
    return 1;
  }
  else if (ffitype == FFI_VOLUME)
  {
    chi_mesh::FieldFunctionInterpolationVolume* new_ffi =
      new chi_mesh::FieldFunctionInterpolationVolume;

    cur_hndlr->ffinterpolation_stack.push_back(new_ffi);
    int index = cur_hndlr->ffinterpolation_stack.size()-1;
    chi_log.Log(LOG_ALLVERBOSE_2)
      << "Created Volume Field Function Interpolation";
    lua_pushnumber(L,index);
    return 1;
  }
  else                                              //Fall back
  {
    chi_log.Log(LOG_ALLERROR)
    << "Invalid FFITypeIndex used in chiFFInterpolationCreate.";
    exit(EXIT_FAILURE);
  }

  return 0;
}