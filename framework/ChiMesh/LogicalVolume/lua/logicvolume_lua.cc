#include"chi_lua.h"

#include "../chi_mesh_logicalvolume.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "chi_runtime.h"

#include "chi_log.h"

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

### Examples
Example usage:
\code
-- Sphere at origin radius 1.0
lv1 = chiLogicalVolumeCreate(SPHERE_ORIGIN, 1.0)

-- Sphere centered at (0.1,0.2,0.3) with radius 1.0
lv2 = chiLogicalVolumeCreate(SPHERE, 0.1, 0.2, 0.3, 1.0)

-- Rectangular parallelepiped
xmin = -1.0; xmax = 1.0
ymin = -2.0; ymax = 1.0
zmin = -1000.0, zmax = 1000.0
lv3 = chiLogicalVolumeCreate(RPP, xmin, xmax, ymin, ymax, zmin, zmax)

-- Right Circular Cylinder
basex = 1.0; basey = 0.0; basez = 0.5
upvecx = 0.0; upvecy = 0.0; upvecz = 1.0
R = 1.0
lv4 = chiLogicalVolumeCreate(RPP, basex, basey, basez, upvecx, upvecy, upvecz, R)

-- Surface mesh
lv_surfmesh = chiSurfaceMeshCreate()
chiSurfaceMeshImportFromOBJFile(lv_surfmesh, "MeshFile3D.obj", false)

lv5 = chiLogicalVolumeCreate(SURFACE, lv_surfmesh)

-- Boolean combination
lv6 = chiLogicalVolumeCreate(BOOLEAN, {{true , lv5},  -- inside logical volume 5
                                       {false, lv1},  -- outside logical volume 1
                                       {false, lv2}}) -- outside logical volume 2
\endcode

\return Handle int Handle to the created logical volume.
\ingroup LuaLogicVolumes
\author Jan*/
int chiLogicalVolumeCreate(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  auto& handler = chi_mesh::GetCurrentHandler();

  const int num_args = lua_gettop(L);
  const int type_index = lua_tonumber(L,1);

  //================================================== Sphere at origin
  if (type_index == SPHERE_ORIGIN)
  {
    if (num_args!=2)
    {
      chi::log.Log0Error() << "Incorrect amount of arguments provided "
                                 "for chiMeshCreateLogicalVolume(SO...";
     chi::Exit(EXIT_FAILURE);
    }
    double r = lua_tonumber(L,2);
    auto log_vol = new chi_mesh::SphereLogicalVolume(r);

    chi::logicvolume_stack.emplace_back(log_vol);
    const size_t index = chi::logicvolume_stack.size()-1;
    lua_pushnumber(L,static_cast<lua_Number>(index));
  }

  //================================================== Sphere at arb loc
  else if (type_index == SPHERE)
  {
    if (num_args!=5)
    {
      chi::log.Log0Error() << "Incorrect amount of arguments provided "
                                 "for chiMeshCreateLogicalVolume(S...";
     chi::Exit(EXIT_FAILURE);
    }
    double x = lua_tonumber(L,2);
    double y = lua_tonumber(L,3);
    double z = lua_tonumber(L,4);
    double r = lua_tonumber(L,5);
    auto log_vol = new chi_mesh::SphereLogicalVolume(x,y,z,r);

    chi::logicvolume_stack.emplace_back(log_vol);
    const size_t index = chi::logicvolume_stack.size()-1;
    lua_pushnumber(L,static_cast<lua_Number>(index));
  }

  //================================================== RPP
  else if (type_index == RPP)
  {
    if (num_args!=7)
    {
      chi::log.Log0Error() << "Incorrect amount of arguments provided "
                                 "for chiMeshCreateLogicalVolume(RPP...";
     chi::Exit(EXIT_FAILURE);
    }
    double xmin = lua_tonumber(L,2);
    double xmax = lua_tonumber(L,3);
    double ymin = lua_tonumber(L,4);
    double ymax = lua_tonumber(L,5);
    double zmin = lua_tonumber(L,6);
    double zmax = lua_tonumber(L,7);
    auto log_vol = new chi_mesh::RPPLogicalVolume(xmin,xmax,ymin,ymax,zmin,zmax);

    chi::logicvolume_stack.emplace_back(log_vol);
    const size_t index = chi::logicvolume_stack.size()-1;
    lua_pushnumber(L,static_cast<lua_Number>(index));
  }

  //================================================== RCC
  else if (type_index == RCC)
  {
    if (num_args!=8)
    {
      chi::log.Log0Error() << "Incorrect amount of arguments provided "
                              "for chiMeshCreateLogicalVolume(RCC...";
      chi::Exit(EXIT_FAILURE);
    }
    double x0 = lua_tonumber(L,2);
    double y0 = lua_tonumber(L,3);
    double z0 = lua_tonumber(L,4);
    double vx = lua_tonumber(L,5);
    double vy = lua_tonumber(L,6);
    double vz = lua_tonumber(L,7);
    double r  = lua_tonumber(L,8);
    auto log_vol = new chi_mesh::RCCLogicalVolume(x0,y0,z0,vx,vy,vz,r);

    chi::logicvolume_stack.emplace_back(log_vol);
    const size_t index = chi::logicvolume_stack.size()-1;
    lua_pushnumber(L,static_cast<lua_Number>(index));

    chi::log.Log0Verbose1()
      << "Created RCC Logical volume with x0,y0,z0,vx,vy,vz,r = "
      << x0 << " " << y0 << " " << z0 << " "
      << vx << " " << vy << " " << vz << " "
      << r;
  }
  else if (type_index == SURFACE)
  {
    if (num_args != 2)
      LuaPostArgAmountError("chiMeshCreateLogicalVolume:SURFACE",2,num_args);

    int surf_mesh_hndle = lua_tonumber(L,2);

    auto surf_mesh_ptr = chi::GetStackItemPtr<chi_mesh::SurfaceMesh>(
      chi::surface_mesh_stack, surf_mesh_hndle, fname);

    auto log_vol = new chi_mesh::SurfaceMeshLogicalVolume(surf_mesh_ptr);

    chi::logicvolume_stack.emplace_back(log_vol);
    const size_t index = chi::logicvolume_stack.size()-1;
    lua_pushnumber(L,static_cast<lua_Number>(index));
  }
  //================================================== BOOLEAN
  else if (type_index == BOOLEAN)
  {
    if (num_args%2 != 0)
    {
      chi::log.Log0Error() << "Incorrect amount of arguments provided "
                                 "for chiMeshCreateLogicalVolume(BOOLEAN..."
                                 " Expected pairs of (bool,volumeHandle)";
     chi::Exit(EXIT_FAILURE);
    }

    auto bool_vol = new chi_mesh::BooleanLogicalVolume;

    int num_pairs = num_args/2;
    for (int p=0; p<num_pairs; p++)
    {
      //==================================== Checking first part of pair
      if (not lua_isboolean(L,2*p))
      {
        chi::log.Log0Error() << "chiMeshCreateLogicalVolume(BOOLEAN..."
                                   " argument " << 2*p << " expected to be "
                                   "Boolean. Found not to be";
       chi::Exit(EXIT_FAILURE);
      }
      //==================================== Checking second part of pair
      if (not lua_isnumber(L,2*p+1))
      {
        chi::log.Log0Error() << "chiMeshCreateLogicalVolume(BOOLEAN..."
                                   " argument " << 2*p+1 << " expected to be "
                                   "number. Found not to be";
       chi::Exit(EXIT_FAILURE);
      }
      if (lua_tonumber(L,2*p+1) >=
          static_cast<lua_Number>(chi::logicvolume_stack.size()))
      {
        chi::log.Log0Error() << "chiMeshCreateLogicalVolume(BOOLEAN..."
                                   " argument " << 2*p+1
                                   << " points to non-existent volume.";
       chi::Exit(EXIT_FAILURE);
      }

      bool logic  = lua_toboolean(L,2*p);
      bool handle = lua_tonumber(L,2*p+1);

      auto p_ref_vol = chi::GetStackItemPtr(chi::logicvolume_stack,
                                            handle, fname);

      typedef std::shared_ptr<chi_mesh::LogicalVolume> LogVolPtr;

      std::pair<bool,LogVolPtr> combo(logic,p_ref_vol);

      bool_vol->parts.push_back(combo);
    }

    chi::logicvolume_stack.emplace_back(bool_vol);
    const size_t index = chi::logicvolume_stack.size()-1;
    lua_pushnumber(L,static_cast<lua_Number>(index));
  }

  //================================================== Unrecognized option
  else
  {
    chi::log.Log0Error() << "Unrecognized volume type used in "
                               "chiLogicalVolumeCreate.";
   chi::Exit(EXIT_FAILURE);
  }


  return 1;
}


//###################################################################
/**Evaluates whether a point is within the logical volume.
\param LVHandle int Handle to the logical volume.
\param Point_x double X-coordinate of the point.
\param Point_y double Y-coordinate of the point.
\param Point_z double Z-coordinate of the point.

\return Sense true if inside the logical volume and false if outside.
*/
int chiLogicalVolumePointSense(lua_State* L)
{
  const std::string fname = "chiLogicalVolumePointSense";
  const int num_args = lua_gettop(L);
  if (num_args != 4)
    LuaPostArgAmountError(fname, 4, num_args);

  LuaCheckNilValue(fname, L, 1);

  const int lv_handle = lua_tointeger(L, 1);

  const auto& lv = chi::GetStackItem<chi_mesh::LogicalVolume>(
    chi::logicvolume_stack, lv_handle, fname);

  const chi_mesh::Vector3 point(lua_tonumber(L,2),
                                lua_tonumber(L,3),
                                lua_tonumber(L,4));

  if (lv.Inside(point)) lua_pushboolean(L, true);
  else                  lua_pushboolean(L, false);

  return 1;
}
