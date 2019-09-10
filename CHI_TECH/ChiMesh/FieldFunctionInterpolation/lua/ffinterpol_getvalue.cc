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
/** Gets the value(s) associated with an interpolation.
 *
\param FFIHandle int Handle to the field function interpolation.

###Note:
Currently only the Volume interpolation supports obtaining a value.

\ingroup LuaFFInterpol
\author Jan*/
int chiFFInterpolationGetValue(lua_State *L)
{
  double value = 0.0;

  int numArgs = lua_gettop(L);
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

  if (typeid(*cur_ffi) != typeid(chi_mesh::FieldFunctionInterpolationVolume))
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiFFInterpolationGetValue is currently only supported for "
      << " VOLUME interpolator types.";
    exit(EXIT_FAILURE);
  }

  chi_mesh::FieldFunctionInterpolationVolume* cur_ffi_volume =
    (chi_mesh::FieldFunctionInterpolationVolume*)cur_ffi;

  value = cur_ffi_volume->op_value;

//  double glob_value = 0.0;
//
//  MPI_Allreduce(&value,&glob_value,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  lua_pushnumber(L,value);
  return 1;
}