#include "../../../ChiLua/chi_lua.h"
#include "../../MeshHandler/chi_meshhandler.h"
#include "../../FieldFunctionInterpolation/Slice/chi_ffinter_slice.h"
#include "../../FieldFunctionInterpolation/Line/chi_ffinter_line.h"
#include "../../FieldFunctionInterpolation/Volume/chi_ffinter_volume.h"

#include "chi_runtime.h"

#include "chi_log.h"
extern ChiLog& chi_log;


#define FFI_FIELD_FUNCTION 0

//#############################################################################
/** Creates a new field function interpolation.
 *
\param FFIHandle int Handle to the field function interpolation.
\param PropertyIndex int Type of property to set.

##_

###PropertyIndex\n
ADD_FIELDFUNCTION     = Add field function to interpolation.\n
SLICE_POINT           = Reference point for SLICE type FFIs.\n
SLICE_NORMAL          = Normal of slice plane.\n
SLICE_TANGENT         = Tangent vector of slice plane.\n
SLICE_BINORM          = Binorm vector of slice plane.\n
LINE_FIRSTPOINT   = Line start point.\n
LINE_SECONDPOINT  = Line end point.\n
LINE_NUMBEROFPOINTS = Number of points to put in the line interpolator.
                          Minimum 2.\n
LINE_CUSTOM_ARRAY = Adds custom array to line interpolator.\n
OPERATION  =  Some interpolations support operation types. See OpTypes.\n
LOGICAL_VOLUME = To be followed by a handle to a logical volume to be
                 used by the interpolator.\n

###OpTypes
Basic operations are volume sum, volume average or volume max. The volume
sum is computed from
\f[
Sum = \sum_k \sum_i u_i \int_V N_i .dV.
\f]
The volume average is computed from
\f[
Avg = \frac{\sum_k \sum_i u_i \int_V N_i .dV}{\sum_k \sum_i \int_V N_i .dV}
\f]
The volume max is simply \f$ max_{k,i}(u_i) \f$.\n\n

OP_SUM\n
For volume interpolations, computes the volume integral.\n
\n
OP_AVG\n
For volume interpolations, computes the volume average.\n
OP_MAX\n
For volume interpolations, computes the volume max.\n
\n
A modified version of these operations are also available. Instead of OP_SUM,
OP_AVG and OP_MAX, the user may supply OP_SUM_LUA, OP_AVG_LUA and OP_MAX_LUA
which then needs to be followed by a string value `LuaFunctionName` of a lua
function of the following form:

\code
function LuaFunctionName(ff_value, mat_id)
 ret_val = 0.0;   --Or some computation
 return ret_val
end
\endcode

This code will be called to return a value \f$ f(u_i) \f$ to be used instead of
the field function \f$ u_i \f$.\n

Example:
\code
xwing=2.0
function IntegrateMaterialVolume(ff_value,mat_id)
    return xwing
end
ffi2 = chiFFInterpolationCreate(VOLUME)
curffi = ffi2
chiFFInterpolationSetProperty(curffi,OPERATION,OP_SUM_LUA,"IntegrateMaterialVolume")
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp)

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
print(chiFFInterpolationGetValue(curffi))
\endcode

The code above will return 2.0 times the volume of cells included in the logical
volume `vol0`.


\return Handle int Handle to the created interpolation.
\ingroup LuaFFInterpol
\author Jan*/
int chiFFInterpolationSetProperty(lua_State *L)
{
  int numArgs = lua_gettop(L);
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  //================================================== Get handle to field function
  const size_t ffihandle = lua_tonumber(L,1);

  auto p_ffi = chi::GetStackItemPtr(chi::field_func_interpolation_stack,
                                    ffihandle, __FUNCTION__);

  //================================================== Process properties
  int property = lua_tonumber(L,2);
  //======================================== Check slice properties
  if ((property >= FFI_PROP_SLICEPOINT) && (property <= FFI_PROP_SLICEBINORM))
  {
    if (typeid(*p_ffi) !=
        typeid(chi_mesh::FieldFunctionInterpolationSlice))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Slice property" << property
        << " used in chiFFInterpolationSetProperty but "
           "FFI is not a slice.";
      chi_log.Log(LOG_0) << typeid(*p_ffi).name();
      chi_log.Log(LOG_0) << typeid(chi_mesh::FieldFunctionInterpolationSlice).name();
      exit(EXIT_FAILURE);
    }

  }

  //======================================== Check Line properties
  if ((property >= FFI_LINE_FIRSTPOINT) && (property <= FFI_LINE_NUMBEROFPOINTS))
  {
    if (typeid(*p_ffi) !=
        typeid(chi_mesh::FieldFunctionInterpolationLine))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Line property" << property
        << " used in chiFFInterpolationSetProperty but "
           "FFI is not a line.";
      chi_log.Log(LOG_0) << typeid(*p_ffi).name();
      chi_log.Log(LOG_0) << typeid(chi_mesh::FieldFunctionInterpolationSlice).name();
      exit(EXIT_FAILURE);
    }

  }

  //========================================= Generic
  if (property == FFI_FIELD_FUNCTION)                        //ADD FF
  {
    int ffhandle = lua_tonumber(L,3);
    std::shared_ptr<chi_physics::FieldFunction> cur_ff = chi::GetStackItemPtr(
      chi::fieldfunc_stack, ffhandle, __FUNCTION__);


    p_ffi->field_functions.push_back(cur_ff);
  }
  else if (property == FFI_PROP_SLICEPOINT)               //REF_POINT
  {
    auto& cur_ffi_slice = (chi_mesh::FieldFunctionInterpolationSlice&)*p_ffi;

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi_slice.point = chi_mesh::Vector3(x, y, z);
  }
  else if (property == FFI_PROP_SLICENORMAL)               //NORMAL
  {
    auto& cur_ffi_slice = (chi_mesh::FieldFunctionInterpolationSlice&)*p_ffi;

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi_slice.normal = chi_mesh::Vector3(x, y, z);
    cur_ffi_slice.normal = cur_ffi_slice.normal/
                            cur_ffi_slice.normal.Norm();
  }
  else if (property == FFI_PROP_SLICETANGENT)               //TANGENT
  {
    auto& cur_ffi_slice = (chi_mesh::FieldFunctionInterpolationSlice&)*p_ffi;

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi_slice.tangent = chi_mesh::Vector3(x, y, z);
    cur_ffi_slice.tangent = cur_ffi_slice.tangent/
                            cur_ffi_slice.tangent.Norm();
  }
  else if (property == FFI_PROP_SLICEBINORM)               //BINORM
  {
    auto& cur_ffi_slice = (chi_mesh::FieldFunctionInterpolationSlice&)*p_ffi;

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi_slice.binorm = chi_mesh::Vector3(x, y, z);
    cur_ffi_slice.binorm = cur_ffi_slice.binorm/
                            cur_ffi_slice.binorm.Norm();
  }
  else if (property == FFI_LINE_FIRSTPOINT)
  {
    if (numArgs!=5)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",5,numArgs);

    auto& cur_ffi_line = (chi_mesh::FieldFunctionInterpolationLine&)*p_ffi;

    cur_ffi_line.pi.x = lua_tonumber(L,3);
    cur_ffi_line.pi.y = lua_tonumber(L,4);
    cur_ffi_line.pi.z = lua_tonumber(L,5);

  }
  else if (property == FFI_LINE_SECONDPOINT)
  {
    if (numArgs!=5)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",5,numArgs);

    auto& cur_ffi_line = (chi_mesh::FieldFunctionInterpolationLine&)*p_ffi;

    cur_ffi_line.pf.x = lua_tonumber(L,3);
    cur_ffi_line.pf.y = lua_tonumber(L,4);
    cur_ffi_line.pf.z = lua_tonumber(L,5);
  }
  else if (property == FFI_LINE_NUMBEROFPOINTS)
  {
    if (numArgs!=3)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",3,numArgs);

    auto& cur_ffi_line = (chi_mesh::FieldFunctionInterpolationLine&)*p_ffi;

    int num_points = lua_tonumber(L,3);

    if (num_points<2)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Line property FFI_LINE_NUMBEROFPOINTS"
        << " used in chiFFInterpolationSetProperty. Number of points must"
        << " be greater than or equal to 2.";
      exit(EXIT_FAILURE);
    }
    cur_ffi_line.number_of_points = num_points;
  }
  else if (property == FFI_LINE_CUSTOM_ARRAY)
  {
    if (numArgs!=3)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",3,numArgs);

    auto& cur_ffi_line = (chi_mesh::FieldFunctionInterpolationLine&)*p_ffi;

    if (not lua_istable(L, 3))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Line property FFI_LINE_CUSTOM_ARRAY"
        << " used in chiFFInterpolationSetProperty. Argument 3 is expected "
           "to be an array.";
      exit(EXIT_FAILURE);
    }

    const size_t table_len = lua_rawlen(L,3);

    std::vector<double> new_array(table_len,0.0);
    for (int k=0; k<table_len; ++k)
    {
      lua_pushnumber(L,k+1);
      lua_gettable(L,3);
      new_array[k] = lua_tonumber(L,-1);
      lua_pop(L,1);
    }

    cur_ffi_line.custom_arrays.push_back(new_array);
  }
  else if (property == FFI_PROP_OPERATION)
  {
    if (numArgs!=3 and numArgs!=4)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",3,numArgs);

    if (typeid(*p_ffi) != typeid(chi_mesh::FieldFunctionInterpolationVolume))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Volume property FFI_PROP_OPERATION"
        << " used in chiFFInterpolationSetProperty can only be used with "
        << "Volume type interpolations.";
      exit(EXIT_FAILURE);
    }

    auto& cur_ffi_volume = (chi_mesh::FieldFunctionInterpolationVolume&)*p_ffi;

    int op_type = lua_tonumber(L,3);

    if (!((op_type>=OP_SUM) && (op_type<=OP_MAX_LUA)))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Volume property FFI_PROP_OPERATION"
        << " used in chiFFInterpolationSetProperty. Unsupported OPERATON."
        << " Supported types are OP_AVG and OP_SUM. " << op_type;
      exit(EXIT_FAILURE);
    }

    if ((op_type >= OP_SUM_LUA) and (op_type <= OP_MAX_LUA))
    {
      if (numArgs != 4)
        LuaPostArgAmountError("chiFFInterpolationSetProperty",4,numArgs);

      const char* func_name = lua_tostring(L,4);
      cur_ffi_volume.op_lua_func = std::string(func_name);
    }

    cur_ffi_volume.op_type = op_type;
  }
  else if (property == FFI_PROP_LOGICAL_VOLUME)
  {
    if (numArgs!=3)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",3,numArgs);

    int logvol_hndle = lua_tonumber(L,3);

    auto p_logical_volume = chi::GetStackItemPtr(chi::logicvolume_stack,
                                                 logvol_hndle,
                                                 __FUNCTION__);

    if (typeid(*p_ffi) != typeid(chi_mesh::FieldFunctionInterpolationVolume))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Volume property FFI_PROP_LOGICAL_VOLUME"
        << " used in chiFFInterpolationSetProperty can only be used with "
        << "Volume type interpolations.";
      exit(EXIT_FAILURE);
    }

    auto& cur_ffi_volume = (chi_mesh::FieldFunctionInterpolationVolume&)*p_ffi;

    cur_ffi_volume.logical_volume = p_logical_volume;
  }
  else                                              //Fall back
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid PropertyIndex used in chiFFInterpolationSetProperty.";
    exit(EXIT_FAILURE);
  }

  return 0;
}
