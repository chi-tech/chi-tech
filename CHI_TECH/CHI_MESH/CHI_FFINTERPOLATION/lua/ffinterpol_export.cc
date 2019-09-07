#include "../../../CHI_LUA/chi_lua.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../../ChiPhysics/chi_physics.h"
#include "../chi_ffinterpolation.h"
#include "../Slice/chi_ffinter_slice.h"
#include "../Line/chi_ffinter_line.h"

#include <chi_log.h>

extern CHI_LOG chi_log;
extern ChiPhysics chi_physics_handler;

//###################################################################
/** Export interpolation to python contour plot.
 *
\param FFIHandle int Handle to the field function interpolation.
\param BaseName char Base name to be used for exported files.

\ingroup LuaFFInterpol
\author Jan*/
int chiFFInterpolationExportPython(lua_State* L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  int num_args = lua_gettop(L);

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

  if (typeid(*cur_ffi) == typeid(chi_mesh::FieldFunctionInterpolationSlice))
  {
    chi_mesh::FieldFunctionInterpolationSlice* cur_ffi_slice =
      (chi_mesh::FieldFunctionInterpolationSlice*)cur_ffi;

    std::string base_name = std::string("ZPFFI") + std::to_string(ffihandle);
    if (num_args==2)
    {
      const char* name = lua_tostring(L,2);
      base_name = std::string(name);
    }

    cur_ffi_slice->ExportPython(base_name);
  }

  if (typeid(*cur_ffi) == typeid(chi_mesh::FieldFunctionInterpolationLine))
  {
    chi_mesh::FieldFunctionInterpolationLine* cur_ffi_line =
      (chi_mesh::FieldFunctionInterpolationLine*)cur_ffi;

    std::string base_name = std::string("ZLFFI") + std::to_string(ffihandle);
    if (num_args==2)
    {
      const char* name = lua_tostring(L,2);
      base_name = std::string(name);
    }

    cur_ffi_line->ExportPython(base_name);
  }


  return 0;
}

