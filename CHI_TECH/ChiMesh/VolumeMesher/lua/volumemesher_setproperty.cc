#include "../../../ChiLua/chi_lua.h"
#include <iostream>
#include "../Predefined2D/volmesher_predefined2d.h"

#include "../../MeshHandler/chi_meshhandler.h"
#include "../../SurfaceMesher/surfacemesher.h"
#include "../../VolumeMesher/Extruder/volmesher_extruder.h"

#include <chi_log.h>
extern ChiLog& chi_log;

//#############################################################################
/** Sets a volume mesher property.

\param PropertyIndex int Index of the property to change. See below
\param PropertyValue varying Value of the property.

##_

###PropertyIndex:
 FORCE_POLYGONS = <B>PropertyValue:[bool]</B> Forces the 2D Meshing to use
                  polygon cells even if the underlying surface mesh is
                  triangles. Expects a boolean value.\n
 MESH_GLOBAL = <B>PropertyValue:[bool]</B> Generate/Read the full mesh at each
               location. Expects a boolean value [Default=true].
 PARTITION_TYPE = <B>PartitionType</B>. See below.
 EXTRUSION_LAYER = <B>PropertyValue:[double,(int),(char)]</B> Adds a layer to the
                   extruder volume mesher if it exists.
                   Expects 1 required parameter, the layer height, followed by 2 optional
                   parameters: number of subdivisions (defaults to 1), and layer id (char)(defaults
                   to nothing). Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```. \n
 CUT_X = Adds a cut at the given x-value. Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 CUT_Y = Adds a cut at the given y-value. Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 CUT_Z = Adds a cut at the given z-value. Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 PARTITION_X   = <B>PropertyValue:[int]</B> Number of partitions in X.
                    Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 PARTITION_Y   = <B>PropertyValue:[int]</B> Number of partitions in Y.
                    Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 PARTITION_Z   = <B>PropertyValue:[int]</B> Number of partitions in Z.
                    Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 MATID_FROMLOGICAL = <B>LogicalVolumeHandle:[int],Mat_id:[int],
                     Sense:[bool](Optional, default:true)</B> Sets the material
                     id of cells that meet the sense requirement for the given
                     logical volume.\n
 BNDRYID_FROMLOGICAL = <B>LogicalVolumeHandle:[int],Bndry_id:[int],
                     Sense:[bool](Optional, default:true)</B> Sets the cell
                     boundary id to the specified value for cells
                     that meet the sense requirement for the given
                     logical volume.\n

## _

### PartitionType
Can be any of the following:
 - KBA_STYLE_XY
 - KBA_STYLE_XYZ
 - PARMETIS

\ingroup LuaVolumeMesher
\author Jan*/
int chiVolumeMesherSetProperty(lua_State *L)
{
  //============================================= Get current mesh handler
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Get property index
  int num_args = lua_gettop(L);
  int property_index = lua_tonumber(L,1);

  typedef chi_mesh::VolumeMesherProperty VMP;

  //============================================= Selects property
  if (property_index == VMP::FORCE_POLYGONS)
  {
    bool force_condition = lua_toboolean(L,2);
    cur_hndlr->volume_mesher->options.force_polygons = force_condition;
  }

  else if (property_index == VMP::MESH_GLOBAL)
  {
    bool mesh_global = lua_toboolean(L,2);
    cur_hndlr->volume_mesher->options.mesh_global = mesh_global;
  }

  else if (property_index == VMP::PARTITION_Z)
  {
    int pz = lua_tonumber(L,2);
    cur_hndlr->volume_mesher->options.partition_z = pz;
    chi_log.Log(LOG_ALLVERBOSE_1)
      << "Partition z set to " << pz;
  }
  else if (property_index == VMP::PARTITION_Y)
  {
    int p = lua_tonumber(L,2);
    cur_hndlr->surface_mesher->partitioning_y = p;
    chi_log.Log(LOG_ALLVERBOSE_1)
      << "Partition y set to " << p;
  }
  else if (property_index == VMP::PARTITION_X)
  {
    int p = lua_tonumber(L,2);
    cur_hndlr->surface_mesher->partitioning_x = p;
    chi_log.Log(LOG_ALLVERBOSE_1)
      << "Partition x set to " << p;
  }
  else if (property_index == VMP::CUTS_Z)
  {
    double p = lua_tonumber(L,2);
    cur_hndlr->volume_mesher->zcuts.push_back(p);
  }
  else if (property_index == VMP::CUTS_Y)
  {
    double p = lua_tonumber(L,2);
    cur_hndlr->surface_mesher->ycuts.push_back(p);
  }
  else if (property_index == VMP::CUTS_X)
  {
    double p = lua_tonumber(L,2);
    cur_hndlr->surface_mesher->xcuts.push_back(p);
  }
  else if (property_index == VMP::PARTITION_TYPE)
  {
    int p = lua_tonumber(L,2);
    if (p >= chi_mesh::VolumeMesher::PartitionType::KBA_STYLE_XY and
        p <= chi_mesh::VolumeMesher::PartitionType::PARMETIS)

    cur_hndlr->volume_mesher->options.partition_type =
      (chi_mesh::VolumeMesher::PartitionType)p;
  }

  else if (property_index == VMP::EXTRUSION_LAYER)
  {
    if (typeid(*cur_hndlr->volume_mesher) == typeid(chi_mesh::VolumeMesherExtruder))
    {
      auto mesher = (chi_mesh::VolumeMesherExtruder*)cur_hndlr->volume_mesher;

      double layer_height = lua_tonumber(L,2);
      int    subdivisions = 1;

      if (num_args>=3)
      {
        subdivisions = lua_tonumber(L,3);
      }
      auto new_layer = new MeshLayer;
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

  else if (property_index == VMP::MATID_FROMLOGICAL)
  {
    if (!((num_args == 3) || (num_args == 4)))
    {
      chi_log.Log(LOG_ALLERROR) << "Invalid amount of arguments used for "
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
      chi_log.Log(LOG_ALLERROR) << "Invalid logical volume specified in "
                                 "chiVolumeMesherSetProperty("
                                 "MATID_FROMLOGICAL...";
      exit(EXIT_FAILURE);
    }

    chi_mesh::LogicalVolume* volume_ptr =
      cur_hndlr->logicvolume_stack[volume_hndl];
    cur_hndlr->volume_mesher->SetMatIDFromLogical(volume_ptr,sense,mat_id);
  }

  else if (property_index == VMP::BNDRYID_FROMLOGICAL)
  {
    if (!((num_args == 3) || (num_args == 4)))
    {
      chi_log.Log(LOG_ALLERROR) << "Invalid amount of arguments used for "
                                 "chiVolumeMesherSetProperty("
                                 "BNDRYID_FROMLOGICAL...";
      exit(EXIT_FAILURE);
    }
    int volume_hndl = lua_tonumber(L,2);
    int bndry_id = lua_tonumber(L,3);
    int sense = true;
    if (num_args==4) sense = lua_toboolean(L,4);

    if (volume_hndl >= cur_hndlr->logicvolume_stack.size())
    {
      chi_log.Log(LOG_ALLERROR) << "Invalid logical volume specified in "
                                 "chiVolumeMesherSetProperty("
                                 "BNDRYID_FROMLOGICAL...";
      exit(EXIT_FAILURE);
    }

    chi_mesh::LogicalVolume* volume_ptr =
      cur_hndlr->logicvolume_stack[volume_hndl];
    cur_hndlr->volume_mesher->SetBndryIDFromLogical(volume_ptr,sense,bndry_id);
  }
  else
  {
    chi_log.Log(LOG_ALLERROR) << "Invalid property specified in call to "
                               "chiVolumeMesherSetProperty().";
    exit(EXIT_FAILURE);
  }

  return 0;
}