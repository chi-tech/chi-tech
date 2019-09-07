#include "../../../CHI_LUA/chi_lua.h"
#include <iostream>
#include "../Linemesh1D/volmesher_linemesh1d.h"
#include "../Predefined2D/volmesher_predefined2d.h"
#include "../Extruder/volmesher_extruder.h"

#include "../../CHI_MESHHANDLER/chi_meshhandler.h"

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
extern ChiLog chi_log;


//#############################################################################
/** Creates a new volume mesher.
 *
\param Type int Volume Remesher type.

Remesher types:\n
 VOLUMEMESHER_LINEMESH1D = Creates 1D slab cells from a linemesh.\n
 VOLUMEMESHER_PREDEFINED2D = No remeshing is performed.\n
 VOLUMEMESHER_EXTRUDER = Extruder the first surface mesh found.\n

\ingroup LuaVolumeMesher
\author Jan*/
int chiVolumeMesherCreate(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  int type = lua_tonumber(L,1);

  chi_mesh::VolumeMesher* new_mesher;
  if (type==VOLUMEMESHER_LINEMESH1D)  //VOLUMEMESHER_PREDEFINED2D
  {
    new_mesher = new chi_mesh::VolumeMesherLinemesh1D;
  }
  else if (type==VOLUMEMESHER_PREDEFINED2D)  //VOLUMEMESHER_PREDEFINED2D
  {
    new_mesher = new chi_mesh::VolumeMesherPredefined2D;
  }
  else if (type==VOLUMEMESHER_EXTRUDER)  //VOLUMEMESHER_EXTRUDER
  {
    new_mesher = new chi_mesh::VolumeMesherExtruder;
  }
  else
  {
    chi_log.Log(LOG_0ERROR) << "Invalid Volume mesher type in function "
                               "chiVolumeMesherCreate. Allowed options are"
                               "VOLUMEMESHER_PREDEFINED2D or "
                               "VOLUMEMESHER_EXTRUDER";
    exit(EXIT_FAILURE);
  }

  cur_hndlr->volume_mesher = new_mesher;

  chi_log.Log(LOG_ALLVERBOSE_2)
    << "chiVolumeMesherCreate: Volume remesher created."
    << std::endl;

  return 0;
}