#include "../../../CHI_LUA/chi_lua.h"
#include <iostream>
#include "../Predefined2D/volmesher_predefined2d.h"

#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../CHI_VOLUMEMESHER/Extruder/volmesher_extruder.h"

#include <chi_log.h>
extern CHI_LOG chi_log;

//#############################################################################
/** Sets a volume mesher property.

\param PropertyIndex int Index of the property to change. See below
\param PropertyValue varying Value of the property.

##_

###PropertyIndex:
 FORCE_POLYGONS = <B>PropertyValue:[bool]</B> Forces the 2D Meshing to use
                  polygon cells even if the
 underlying surface mesh is triangles. Expects a boolean value.\n
 EXTRUSION_LAYER = <B>PropertyValue:[double,(int),(char)]</B> Adds a layer to the
                   extruder volume mesher if it exists.
                   Expects 1 required parameter, the layer height, followed by 2 optional
                   parameters: number of subdivisions (defaults to 1), and layer id (char)(defaults
                   to nothing).\n
 MESH_GLOBAL = <B>PropertyValue:[bool]</B> Generate/Read the full mesh at each
               location. Expects a boolean value [Default=true].\n
 PARTITION_Z = <B>PropertyValue:[int]</B> Sets the pz partitioning parameters
                 for the z direction.\n
 MATID_FROMLOGICAL = <B>LogicalVolumeHandle:[int],Mat_id:[int],
                     Sense:[bool](Optional, default:true)</B> Sets the material
                     id of cells that meet the sense requirement for the given
                     logical volume.


\ingroup LuaVolumeMesher
\author Jan*/
int chiVolumeMesherSetProperty(lua_State *L)
{
  //============================================= Get current mesh handler
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Get property index
  int num_args = lua_gettop(L);
  int property_index = lua_tonumber(L,1);

  //============================================= Selects property
  if (property_index == 1) //FORCE_POLYGONS
  {
    bool force_condition = lua_toboolean(L,2);
    cur_hndlr->volume_mesher->options.force_polygons = force_condition;
  }

  if (property_index == 2) //MESH_GLOBAL
  {
    bool mesh_global = lua_toboolean(L,2);
    cur_hndlr->volume_mesher->options.mesh_global = mesh_global;
  }

  if (property_index == 3) //PARTITION_Z
  {
    int pz = lua_tonumber(L,2);
    cur_hndlr->volume_mesher->options.partition_z = pz;
    chi_log.Log(LOG_ALLVERBOSE_2)
      << "Partition z set to " << pz;
  }


  if (property_index == 10) //EXTRUSION_LAYER
  {
    if (typeid(*cur_hndlr->volume_mesher) == typeid(chi_mesh::VolumeMesherExtruder))
    {
      chi_mesh::VolumeMesherExtruder* mesher =
        (chi_mesh::VolumeMesherExtruder*)cur_hndlr->volume_mesher;

      double layer_height = lua_tonumber(L,2);
      int    subdivisions = 1;

      if (num_args>=3)
      {
        subdivisions = lua_tonumber(L,3);
      }
      MESH_LAYER* new_layer = new MESH_LAYER;
      new_layer->height = layer_height;
      new_layer->sub_divisions = subdivisions;

      if (num_args==4)
      {
        new_layer->name = std::string(lua_tostring(L,4));
      }
      mesher->input_layers.push_back(new_layer);
    }
    else
    {
      fprintf(stdout,"VolumeMesherExtruder is not the current volume mesher"
                     " therefore the z-layer property is ignored.\n");
    }

  }

  if (property_index == 11) //MATID_FROMLOGICAL
  {
    if (!((num_args == 3) || (num_args == 4)))
    {
      chi_log.Log(LOG_0ERROR) << "Invalid amount of arguments used for "
                                 "chiVolumeMesherSetProperty("
                                 "MATID_FROMLOGICAL...";
      exit(EXIT_FAILURE);
    }
    int volume_hndl = lua_tonumber(L,2);
    int mat_id = lua_tonumber(L,3);
    int sense = true;
    if (num_args==4) sense = lua_toboolean(L,4);

    if (volume_hndl >= cur_hndlr->logicvolume_stack.size())
    {
      chi_log.Log(LOG_0ERROR) << "Invalid logical volume specified in "
                                 "chiVolumeMesherSetProperty("
                                 "MATID_FROMLOGICAL...";
      exit(EXIT_FAILURE);
    }

    chi_mesh::LogicalVolume* volume_ptr =
      cur_hndlr->logicvolume_stack[volume_hndl];
    cur_hndlr->volume_mesher->SetMatIDFromLogical(volume_ptr,sense,mat_id);
  }

  return 0;
}