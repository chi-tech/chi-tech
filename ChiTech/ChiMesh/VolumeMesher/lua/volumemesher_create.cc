#include "ChiLua/chi_lua.h"
#include "ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h"
#include "ChiMesh/VolumeMesher/PredefinedUnpartitioned/volmesher_predefunpart.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

#include <iostream>

#include "chi_runtime.h"

#include "chi_log.h"
extern ChiLog& chi_log;


//#############################################################################
/** Creates a new volume mesher.
 *
\param Type int Volume Remesher type.
\param OtherArgs varying Additional arguments depending on `Type`.

Remesher types:\n
 VOLUMEMESHER_EXTRUDER = Creates an extruded mesh from a 2D template mesh.
 Requires two additional arguments. `TemplateType` and `handle`. See below.\n
 VOLUMEMESHER_UNPARTITIONED = Create the mesh from the latest UnpartitionedMesh.
 Requires a single additional argument, `handle`, which is a handle to
 a valid unpartitioned mesh.\n

##_

###Extruder parameters

When the mesher type is specified to be VOLUMEMESHER_EXTRUDER then two
additional arguments are required. `TemplateType` and `handle`.\n

- `TemplateType` can be either `ExtruderTemplateType.SURFACE_MESH` or
 `ExtruderTemplateType.UNPARTITIONED_MESH`.\n
- `handle` is a handle to the template mesh. When `TemplateType` is
 set to `ExtruderTemplateType.SURFACE_MESH` then this must be a handle to a valid
 surface mesh. Similarly, when `TemplateType` is set to
 `ExtruderTemplateType.UNPARTITIONED_MESH` then the handle must point to a valid
 unpartitioned mesh.

\ingroup LuaVolumeMesher
\author Jan*/
int chiVolumeMesherCreate(lua_State *L)
{
  const std::string fname = __FUNCTION__;

  //============================================= Arguments check
  const int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  //============================================= Mesher type
  const int mesher_type = lua_tonumber(L, 1);

  chi_mesh::VolumeMesher* new_mesher;

  if (mesher_type == chi_mesh::VolumeMesherType::EXTRUDER)
  {
    if (num_args != 3)
    {
      chi_log.Log(LOG_ALLERROR)
        << fname + ": "
           "When specifying VOLUMEMESHER_EXTRUDER, the template type and "
           "handle must also be supplied.";
      exit(EXIT_FAILURE);
    }

    LuaCheckNilValue(fname,L,2);
    LuaCheckNilValue(fname,L,3);

    int template_type   = lua_tonumber(L,2);
    int template_handle = lua_tonumber(L,3);

    const auto SURFACE_MESH_TEMPLATE =
      chi_mesh::VolumeMesherExtruder::TemplateType::SURFACE_MESH;
    const auto UNPART_MESH_TEMPLATE =
      chi_mesh::VolumeMesherExtruder::TemplateType::UNPARTITIONED_MESH;

    auto& handler = chi_mesh::GetCurrentHandler();

    if      (template_type == (int)SURFACE_MESH_TEMPLATE)
    {
      auto surface_mesh_ptr = chi::GetStackItemPtr(chi::surface_mesh_stack,
                                                   template_handle, fname);

      new_mesher = new chi_mesh::VolumeMesherExtruder(surface_mesh_ptr);
    }
    else if (template_type == (int)UNPART_MESH_TEMPLATE)
    {
      auto p_umesh = chi::GetStackItemPtr(chi::unpartitionedmesh_stack,
                                          template_handle, fname);

      new_mesher = new chi_mesh::VolumeMesherExtruder(p_umesh);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "In call to " << __FUNCTION__ << ". Invalid template type specified.";
      exit(EXIT_FAILURE);
    }


  }
  else if (mesher_type == chi_mesh::VolumeMesherType::UNPARTITIONED)
  {
    if (num_args != 2)
    {
      chi_log.Log(LOG_ALLERROR)
        << fname + ": "
                   "When specifying VOLUMEMESHER_UNPARTITIONED, the "
                   "handle must also be supplied.";
      exit(EXIT_FAILURE);
    }

    LuaCheckNilValue(fname,L,2);
    const int template_handle = lua_tonumber(L,2);

    auto p_umesh = chi::GetStackItemPtr(chi::unpartitionedmesh_stack,
                                        template_handle, fname);

    new_mesher = new chi_mesh::VolumeMesherPredefinedUnpartitioned(p_umesh);
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