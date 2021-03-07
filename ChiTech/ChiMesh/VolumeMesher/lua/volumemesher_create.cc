#include "../../../ChiLua/chi_lua.h"
#include "../Predefined2D/volmesher_predefined2d.h"
#include "../Extruder/volmesher_extruder.h"
#include "../PredefinedUnpartitioned/volmesher_predefunpart.h"

#include "../../MeshHandler/chi_meshhandler.h"

#include <iostream>

/** \defgroup LuaVolumeMesher Volume Meshers
 * \ingroup LuaMesh
 *
 * ## General Concepts
 *
 * A volume mesher generates cells (chi_mesh::Cell) which can be
 * either 2D polygons (chi_mesh::CellPolygon) or 3D polyhedrons
 * (chi_mesh::CellPolyhedron). All cell objects are pushed into
 * a vector located in a chi_mesh::MeshContinuum after creation.
 * Right now only a single continuum can be operated on and this is normally
 * by means of attaching a region to a solver. \n
 *
 * \n
 * Right now 2D meshes can be partitioned but the full 2D cell geometry
 * exists on each processor. This is mostly because of the way the extruder
 * works. This means that for a 2D problem a continuum will have ALL of the
 * cells (as fully defined cells), in every process, in the member
 * vector "cells" of the continuum.
 * The same can not be said for 3D extruded meshes where each process has all
 * the nodes but if cells are not local then only placeholders are uploaded.
 * Placeholders are basically the base class chi_mesh::Cell and contains only
 * the cell's partition_id and its centroid. The concept that these placeholders
 * are uploaded allows a sweeping order to figure out dependencies. A similar
 * strategy would have to be devised for using third-party meshes.\n
 *
 * \n
 * Global mesh references are maintained in chi_mesh::MeshContinuum::cells.
 * This contains the actual mesh object. Local indices are stored in
 * chi_mesh::MeshContinuum::local_cell_glob_indices and are the global indices
 * of local cells. Conversely the local indices, given a global index,
 * are stored in chi_mesh::MeshContinuum::glob_cell_local_indices.
 *
 * ## Extruder Mesher
*/

#include <chi_log.h>
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

  if (type==chi_mesh::VolumeMesherType::PREDEFINED2D)
  {
    new_mesher = new chi_mesh::VolumeMesherPredefined2D;
  }
  else if (type==chi_mesh::VolumeMesherType::EXTRUDER)
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

    auto handler = chi_mesh::GetCurrentHandler();

    if      (template_type == (int)SURFACE_MESH_TEMPLATE)
    {
      chi_mesh::SurfaceMesh* surfaceMesh;
      try {
        surfaceMesh = handler->surface_mesh_stack.at(template_handle);
      }
      catch (const std::out_of_range& o)
      {
        chi_log.Log(LOG_ALLERROR)
          << "In call to " << __FUNCTION__ << ", with template type"
          << " ExtruderTemplateType.SURFACE_MESH. Invalid handle "
          << template_handle << " to surface mesh.";
        exit(EXIT_FAILURE);
      }
      new_mesher = new chi_mesh::VolumeMesherExtruder(surfaceMesh);
    }
    else if (template_type == (int)UNPART_MESH_TEMPLATE)
    {
      chi_mesh::UnpartitionedMesh* unpartitionedMesh;
      try {
        unpartitionedMesh = handler->unpartitionedmesh_stack.at(template_handle);
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
                               "VOLUMEMESHER_PREDEFINED2D, "
                               "VOLUMEMESHER_EXTRUDER or "
                               "VOLUMEMESHER_UNPARTITIONED";
    exit(EXIT_FAILURE);
  }

  auto cur_hndlr = chi_mesh::GetCurrentHandler();
  cur_hndlr->volume_mesher = new_mesher;

  chi_log.Log(LOG_ALLVERBOSE_2)
    << "chiVolumeMesherCreate: Volume mesher created."
    << std::endl;

  return 0;
}