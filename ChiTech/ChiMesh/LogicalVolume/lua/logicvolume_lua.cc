#include"chi_lua.h"
#include <iostream>

#include "../chi_mesh_logicalvolume.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include <chi_log.h>

extern ChiLog& chi_log;

/** \defgroup LuaLogicVolumes Logical Volumes
 * \ingroup LuaMesh*/

//###################################################################
/** Creates a logical volume.

\param TypeIndex int Volume type.
\param Values varying Parameters.

##_

###TypeIndex\n
SPHERE_ORIGIN =  Sphere at the origin. [Requires: R]\n
SPHERE  =  Sphere at user supplied location. [Requires: x,y,z,R]\n
RPP=  Rectangular ParalleliPiped. [Requires: xmin,xmax,ymin,ymax,zmin,zmax]\n
RCC=  Right Circular Cylinder. [Requires: x0,y0,z0, vx,vy,vz and R]\n
SURFACE= Logical volume determined from a surface mesh. [Requires: a handle
    to a surface mesh]\n
BOOLEAN= Boolean combination of other volumes.
 [Requires pairs of values: bool,volumeHandle]\n



\return Handle int Handle to the created logical volume.
\ingroup LuaLogicVolumes
\author Jan*/
int chiLogicalVolumeCreate(lua_State *L)
{
  auto& handler = chi_mesh::GetCurrentHandler();


  int num_args = lua_gettop(L);
  int type_index = lua_tonumber(L,1);

  //================================================== Sphere at origin
  if (type_index == SPHERE_ORIGIN)
  {
    if (num_args!=2)
    {
      chi_log.Log(LOG_0ERROR) << "Incorrect amount of arguments provided "
                                 "for chiMeshCreateLogicalVolume(SO...";
      exit(EXIT_FAILURE);
    }
    double r = lua_tonumber(L,2);
    chi_mesh::SphereLogicalVolume* log_vol =
      new chi_mesh::SphereLogicalVolume(r);

    handler.logicvolume_stack.push_back(log_vol);
    lua_pushnumber(L,handler.logicvolume_stack.size()-1);
  }

  //================================================== Sphere at arb loc
  else if (type_index == SPHERE)
  {
    if (num_args!=5)
    {
      chi_log.Log(LOG_0ERROR) << "Incorrect amount of arguments provided "
                                 "for chiMeshCreateLogicalVolume(S...";
      exit(EXIT_FAILURE);
    }
    double x = lua_tonumber(L,2);
    double y = lua_tonumber(L,3);
    double z = lua_tonumber(L,4);
    double r = lua_tonumber(L,5);
    chi_mesh::SphereLogicalVolume* log_vol =
      new chi_mesh::SphereLogicalVolume(x,y,z,r);

    handler.logicvolume_stack.push_back(log_vol);
    lua_pushnumber(L,handler.logicvolume_stack.size()-1);
  }

  //================================================== RPP
  else if (type_index == RPP)
  {
    if (num_args!=7)
    {
      chi_log.Log(LOG_0ERROR) << "Incorrect amount of arguments provided "
                                 "for chiMeshCreateLogicalVolume(RPP...";
      exit(EXIT_FAILURE);
    }
    double xmin = lua_tonumber(L,2);
    double xmax = lua_tonumber(L,3);
    double ymin = lua_tonumber(L,4);
    double ymax = lua_tonumber(L,5);
    double zmin = lua_tonumber(L,6);
    double zmax = lua_tonumber(L,7);
    chi_mesh::RPPLogicalVolume* log_vol =
      new chi_mesh::RPPLogicalVolume(xmin,xmax,ymin,ymax,zmin,zmax);

    handler.logicvolume_stack.push_back(log_vol);
    lua_pushnumber(L,handler.logicvolume_stack.size()-1);
  }

  //================================================== RCC
  else if (type_index == RCC)
  {
    if (num_args!=8)
    {
      chi_log.Log(LOG_0ERROR) << "Incorrect amount of arguments provided "
                                 "for chiMeshCreateLogicalVolume(RCC...";
      exit(EXIT_FAILURE);
    }
    double x0 = lua_tonumber(L,2);
    double y0 = lua_tonumber(L,3);
    double z0 = lua_tonumber(L,4);
    double vx = lua_tonumber(L,5);
    double vy = lua_tonumber(L,6);
    double vz = lua_tonumber(L,7);
    double r    = lua_tonumber(L,7);
    chi_mesh::RCCLogicalVolume* log_vol =
      new chi_mesh::RCCLogicalVolume(x0,y0,z0,vx,vy,vz,r);

    handler.logicvolume_stack.push_back(log_vol);
    lua_pushnumber(L,handler.logicvolume_stack.size()-1);

    chi_mesh::Vector3 point(-0.5, 0.0, 0.1);
    printf("MATRIX %d\n", log_vol->Inside(point));
  }
  else if (type_index == SURFACE)
  {
    if (num_args != 2)
      LuaPostArgAmountError("chiMeshCreateLogicalVolume:SURFACE",2,num_args);

    int surf_mesh_hndle = lua_tonumber(L,2);

    chi_mesh::SurfaceMesh* surf_mesh;
    try {
      surf_mesh = handler.surface_mesh_stack.at(surf_mesh_hndle);
    }
    catch(const std::out_of_range& o)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Invalid handle to surface mesh specified in call to "
           "chiMeshCreateLogicalVolume:SURFACE";
      exit(EXIT_FAILURE);
    }

    chi_mesh::SurfaceMeshLogicalVolume* surf_vol =
      new chi_mesh::SurfaceMeshLogicalVolume(surf_mesh);

    handler.logicvolume_stack.push_back(surf_vol);
    lua_pushnumber(L,handler.logicvolume_stack.size()-1);
  }
  //================================================== BOOLEAN
  else if (type_index == BOOLEAN)
  {
    if (num_args%2 != 0)
    {
      chi_log.Log(LOG_0ERROR) << "Incorrect amount of arguments provided "
                                 "for chiMeshCreateLogicalVolume(BOOLEAN..."
                                 " Expected pairs of (bool,volumeHandle)";
      exit(EXIT_FAILURE);
    }

    chi_mesh::BooleanLogicalVolume* bool_vol =
      new chi_mesh::BooleanLogicalVolume;

    int num_pairs = num_args/2;
    for (int p=0; p<num_pairs; p++)
    {
      //==================================== Checking first part of pair
      if (not lua_isboolean(L,2*p))
      {
        chi_log.Log(LOG_0ERROR) << "chiMeshCreateLogicalVolume(BOOLEAN..."
                                   " argument " << 2*p << " expected to be "
                                   "Boolean. Found not to be";
        exit(EXIT_FAILURE);
      }
      //==================================== Checking second part of pair
      if (not lua_isnumber(L,2*p+1))
      {
        chi_log.Log(LOG_0ERROR) << "chiMeshCreateLogicalVolume(BOOLEAN..."
                                   " argument " << 2*p+1 << " expected to be "
                                   "number. Found not to be";
        exit(EXIT_FAILURE);
      }
      if (lua_tonumber(L,2*p+1)>=handler.logicvolume_stack.size())
      {
        chi_log.Log(LOG_0ERROR) << "chiMeshCreateLogicalVolume(BOOLEAN..."
                                   " argument " << 2*p+1
                                   << " points to non-existent volume.";
        exit(EXIT_FAILURE);
      }

      bool logic  = lua_toboolean(L,2*p);
      bool handle = lua_tonumber(L,2*p+1);

      chi_mesh::LogicalVolume* ref_vol = handler.logicvolume_stack[handle];
      std::pair<bool,chi_mesh::LogicalVolume*> combo(logic,ref_vol);

      bool_vol->parts.push_back(combo);
    }

    handler.logicvolume_stack.push_back(bool_vol);
  }

  //================================================== Unrecognized option
  else
  {
    chi_log.Log(LOG_0ERROR) << "Unrecognized volume type used in "
                               "chiLogicalVolumeCreate.";
    exit(EXIT_FAILURE);
  }


  return 1;
}
