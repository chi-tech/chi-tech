#include "../../../ChiLua/chi_lua.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../CHI_FFINTERPOLATION/Slice/chi_ffinter_slice.h"
#include "../../CHI_FFINTERPOLATION/Line/chi_ffinter_line.h"
#include "../../CHI_FFINTERPOLATION/Volume/chi_ffinter_volume.h"
#include <ChiMesh/CHI_LOGICALVOLUME/chi_mesh_logicalvolume.h>
#include "../../../ChiPhysics/chi_physics.h"

#include <chi_log.h>

extern ChiLog chi_log;
extern ChiPhysics chi_physics_handler;

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
OPERATION  =  Some interpolations support operation types. See OpTypes.\n
LOGICAL_VOLUME = To be followed by a handle to a logical volume to be
                 used by the interpolator.\n

###OpTypes
OP_SUM\n
For volume interpolations, computes the volume integral.\n
\n
OP_AVG\n
For volume interpolations, computes the volume average.\n
OP_MAX\n
For volume interpolations, computes the volume max.\n



\return Handle int Handle to the created interpolation.
\ingroup LuaFFInterpol
\author Jan*/
int chiFFInterpolationSetProperty(lua_State *L)
{
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

  //================================================== Process properties
  int property = lua_tonumber(L,2);
  //======================================== Check slice properties
  if ((property >= FFI_PROP_SLICEPOINT) && (property <= FFI_PROP_SLICEBINORM))
  {
    if (typeid(*cur_ffi) !=
        typeid(chi_mesh::FieldFunctionInterpolationSlice))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Slice property" << property
        << " used in chiFFInterpolationSetProperty but "
           "FFI is not a slice.";
      chi_log.Log(LOG_0) << typeid(*cur_ffi).name();
      chi_log.Log(LOG_0) << typeid(chi_mesh::FieldFunctionInterpolationSlice).name();
      exit(EXIT_FAILURE);
    }

  }

  //======================================== Check Line properties
  if ((property >= FFI_LINE_FIRSTPOINT) && (property <= FFI_LINE_NUMBEROFPOINTS))
  {
    if (typeid(*cur_ffi) !=
        typeid(chi_mesh::FieldFunctionInterpolationLine))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Line property" << property
        << " used in chiFFInterpolationSetProperty but "
           "FFI is not a line.";
      chi_log.Log(LOG_0) << typeid(*cur_ffi).name();
      chi_log.Log(LOG_0) << typeid(chi_mesh::FieldFunctionInterpolationSlice).name();
      exit(EXIT_FAILURE);
    }

  }

  //========================================= Generic
  if (property == FFI_FIELD_FUNCTION)                        //ADD FF
  {
    int ffhandle = lua_tonumber(L,3);
    chi_physics::FieldFunction* cur_ff;
    try {
      cur_ff = chi_physics_handler.fieldfunc_stack.at(ffhandle);
    }
    catch(std::out_of_range o)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Invalid field function handle in chiFFInterpolationSetProperty.";
      exit(EXIT_FAILURE);
    }

    cur_ffi->field_functions.push_back(cur_ff);
  }
  else if (property == FFI_PROP_SLICEPOINT)               //REF_POINT
  {
    chi_mesh::FieldFunctionInterpolationSlice* cur_ffi_slice =
      (chi_mesh::FieldFunctionInterpolationSlice*)cur_ffi;

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi_slice->point = chi_mesh::Vector(x,y,z);
  }
  else if (property == FFI_PROP_SLICENORMAL)               //NORMAL
  {
    chi_mesh::FieldFunctionInterpolationSlice* cur_ffi_slice =
      (chi_mesh::FieldFunctionInterpolationSlice*)cur_ffi;

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi_slice->normal = chi_mesh::Vector(x,y,z);
    cur_ffi_slice->normal = cur_ffi_slice->normal/
                            cur_ffi_slice->normal.Norm();
  }
  else if (property == FFI_PROP_SLICETANGENT)               //TANGENT
  {
    chi_mesh::FieldFunctionInterpolationSlice* cur_ffi_slice =
      (chi_mesh::FieldFunctionInterpolationSlice*)cur_ffi;

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi_slice->tangent = chi_mesh::Vector(x,y,z);
    cur_ffi_slice->tangent = cur_ffi_slice->tangent/
                            cur_ffi_slice->tangent.Norm();
  }
  else if (property == FFI_PROP_SLICEBINORM)               //BINORM
  {
    chi_mesh::FieldFunctionInterpolationSlice* cur_ffi_slice =
      (chi_mesh::FieldFunctionInterpolationSlice*)cur_ffi;

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi_slice->binorm = chi_mesh::Vector(x,y,z);
    cur_ffi_slice->binorm = cur_ffi_slice->binorm/
                            cur_ffi_slice->binorm.Norm();
  }
  else if (property == FFI_LINE_FIRSTPOINT)
  {
    if (numArgs!=5)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",5,numArgs);

    chi_mesh::FieldFunctionInterpolationLine* cur_ffi_line =
      (chi_mesh::FieldFunctionInterpolationLine*)cur_ffi;

    cur_ffi_line->pi.x = lua_tonumber(L,3);
    cur_ffi_line->pi.y = lua_tonumber(L,4);
    cur_ffi_line->pi.z = lua_tonumber(L,5);

  }
  else if (property == FFI_LINE_SECONDPOINT)
  {
    if (numArgs!=5)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",5,numArgs);

    chi_mesh::FieldFunctionInterpolationLine* cur_ffi_line =
      (chi_mesh::FieldFunctionInterpolationLine*)cur_ffi;

    cur_ffi_line->pf.x = lua_tonumber(L,3);
    cur_ffi_line->pf.y = lua_tonumber(L,4);
    cur_ffi_line->pf.z = lua_tonumber(L,5);
  }
  else if (property == FFI_LINE_NUMBEROFPOINTS)
  {
    if (numArgs!=3)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",3,numArgs);

    chi_mesh::FieldFunctionInterpolationLine* cur_ffi_line =
      (chi_mesh::FieldFunctionInterpolationLine*)cur_ffi;

    int num_points = lua_tonumber(L,3);

    if (num_points<2)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Line property FFI_LINE_NUMBEROFPOINTS"
        << " used in chiFFInterpolationSetProperty. Number of points must"
        << " be greater than or equal to 2.";
      exit(EXIT_FAILURE);
    }
    cur_ffi_line->number_of_points = num_points;
  }
  else if (property == FFI_PROP_OPERATION)
  {
    if (numArgs!=3)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",3,numArgs);

    if (typeid(*cur_ffi) != typeid(chi_mesh::FieldFunctionInterpolationVolume))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Volume property FFI_PROP_OPERATION"
        << " used in chiFFInterpolationSetProperty can only be used with "
        << "Volume type interpolations.";
      exit(EXIT_FAILURE);
    }

    chi_mesh::FieldFunctionInterpolationVolume* cur_ffi_volume =
      (chi_mesh::FieldFunctionInterpolationVolume*)cur_ffi;

    int op_type = lua_tonumber(L,3);

    if (!((op_type>=OP_SUM) && (op_type<=OP_MAX)))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Volume property FFI_PROP_OPERATION"
        << " used in chiFFInterpolationSetProperty. Unsupported OPERATON."
        << " Supported types are OP_AVG and OP_SUM. " << op_type;
      exit(EXIT_FAILURE);
    }
    cur_ffi_volume->op_type = op_type;
  }
  else if (property == FFI_PROP_LOGICAL_VOLUME)
  {
    if (numArgs!=3)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",3,numArgs);

    int logvol_hndle = lua_tonumber(L,3);

    chi_mesh::LogicalVolume* logvol = nullptr;

    try {
      logvol = cur_hndlr->logicvolume_stack.at(logvol_hndle);
    }
    catch(std::out_of_range o)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Invalid logical volume handle in chiFFInterpolationSetProperty.";
      exit(EXIT_FAILURE);
    }

    if (typeid(*cur_ffi) != typeid(chi_mesh::FieldFunctionInterpolationVolume))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Volume property FFI_PROP_LOGICAL_VOLUME"
        << " used in chiFFInterpolationSetProperty can only be used with "
        << "Volume type interpolations.";
      exit(EXIT_FAILURE);
    }

    chi_mesh::FieldFunctionInterpolationVolume* cur_ffi_volume =
      (chi_mesh::FieldFunctionInterpolationVolume*)cur_ffi;

    cur_ffi_volume->logical_volume = logvol;
  }
  else                                              //Fall back
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid PropertyIndex used in chiFFInterpolationSetProperty.";
    exit(EXIT_FAILURE);
  }

  return 0;
}