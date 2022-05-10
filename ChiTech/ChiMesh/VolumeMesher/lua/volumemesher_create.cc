#include "ChiLua/chi_lua.h"
#include "ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h"
#include "ChiMesh/VolumeMesher/PredefinedUnpartitioned/volmesher_predefunpart.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include <iostream>

#include "chi_runtime.h"

#include "chi_log.h"
extern ChiLog& chi_log;


//#############################################################################
/** Creates a new volume mesher.
 *
\param Type int Volume Remesher type.

Remesher types:\n
 VOLUMEMESHER_PREDEFINED2D = No remeshing is performed.\n
 VOLUMEMESHER_EXTRUDER = Creates an extruded mesh from a 2D template mesh.
 Requires two additional arguments. <I>TemplateType</I> and <I>handle</I>. See below.\n
 VOLUMEMESHER_UNPARTITIONED = Create the mesh from the latest UnpartitionedMesh.\n

##_

###Extruder parameters

When the mesher type is specified to be VOLUMEMESHER_EXTRUDER then two
additional arguments are required. <I>TemplateType</I> and <I>handle</I>.\n

- <I>TemplateType</I> can be either ExtruderTemplateType.SURFACE_MESH or
 ExtruderTemplateType.UNPARTITIONED_MESH.\n
- <I>handle</I> is a handle to the template mesh. When <I>TemplateType</I> is
 set to ExtruderTemplateType.SURFACE_MESH then this must be a handle to a valid
 surface mesh. Similarly, when <I>TemplateType</I> is set to
 ExtruderTemplateType.UNPARTITIONED_MESH then the handle must point to a valid
 unpartitioned mesh.

\ingroup LuaVolumeMesher
\author Jan*/
int chiVolumeMesherCreate(lua_State *L)
{
  //============================================= Arguments check
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  //============================================= Process type
  int type = lua_tonumber(L,1);

  chi_mesh::VolumeMesher* new_mesher;

  if (type==chi_mesh::VolumeMesherType::EXTRUDER)
  {
    if (num_args != 3)
    {
      chi_log.Log(LOG_ALLERROR)
        << "When specifying VOLUMEMESHER_EXTRUDER, the template type and "
           "handle must also be supplied.";
      exit(EXIT_FAILURE);
    }

    LuaCheckNilValue(__FUNCTION__,L,2);
    LuaCheckNilValue(__FUNCTION__,L,3);

    int template_type   = lua_tonumber(L,2);
    int template_handle = lua_tonumber(L,3);

    const auto SURFACE_MESH_TEMPLATE =
      chi_mesh::VolumeMesherExtruder::TemplateType::SURFACE_MESH;
    const auto UNPART_MESH_TEMPLATE =
      chi_mesh::VolumeMesherExtruder::TemplateType::UNPARTITIONED_MESH;

    auto& handler = chi_mesh::GetCurrentHandler();

    if      (template_type == (int)SURFACE_MESH_TEMPLATE)
    {
      auto surface_mesh_ptr = chi::GetStackItemPtr(
        chi::surface_mesh_stack, template_handle, __FUNCTION__);
      new_mesher = new chi_mesh::VolumeMesherExtruder(surface_mesh_ptr);
    }
    else if (template_type == (int)UNPART_MESH_TEMPLATE)
    {
      chi_mesh::UnpartitionedMesh* unpartitionedMesh;
      try {
        unpartitionedMesh = handler.unpartitionedmesh_stack.at(template_handle);
      }
      catch (const std::out_of_range& o)
      {
        chi_log.Log(LOG_ALLERROR)
          << "In call to " << __FUNCTION__ << ", with template type"
          << " ExtruderTemplateType.UNPARTITIONED_MESH. Invalid handle "
          << template_handle << " to unpartitioned mesh.";
        exit(EXIT_FAILURE);
      }
      new_mesher = new chi_mesh::VolumeMesherExtruder(unpartitionedMesh);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "In call to " << __FUNCTION__ << ". Invalid template type specified.";
      exit(EXIT_FAILURE);
    }


  }
  else if (type==chi_mesh::VolumeMesherType::UNPARTITIONED)
  {
    new_mesher = new chi_mesh::VolumeMesherPredefinedUnpartitioned;
  }
  else
  {
    chi_log.Log(LOG_0ERROR) << "Invalid Volume mesher type in function "
                               "chiVolumeMesherCreate. Allowed options are"
                               "VOLUMEMESHER_EXTRUDER or "
                               "VOLUMEMESHER_UNPARTITIONED";
    exit(EXIT_FAILURE);
  }

  auto& cur_hndlr = chi_mesh::GetCurrentHandler();
  cur_hndlr.volume_mesher = new_mesher;

  chi_log.Log(LOG_ALLVERBOSE_2)
    << "chiVolumeMesherCreate: Volume mesher created."
    << std::endl;

  return 0;
}