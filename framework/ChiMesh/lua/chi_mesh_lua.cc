#include "chi_mesh_lua.h"

#include "ChiMesh/OrthoGrids/lua/lua_mesh_orthomacros.h"
#include "ChiMesh/UnpartitionedMesh/lua/unpartition_mesh_lua_utils.h"
#include "ChiMesh/LogicalVolume/lua/logicvolume_lua.h"
#include "ChiMesh/MeshHandler/lua/meshhandler_lua.h"
#include "ChiMesh/SurfaceMesher/lua/surfmesher_lua.h"
#include "ChiMesh/VolumeMesher/lua/volumemesher_lua.h"
#include "ChiMesh/DomainDecomposition/lua/domaindecomp_lua.h"
#include "ChiMesh/MeshCutting/lua/meshcutting_lua.h"
#include "ChiMesh/FieldFunctionInterpolation/lua/ffinterpol_lua.h"
#include "ChiMesh/SurfaceMesh/lua/lua_surface_mesh.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)
#define LUA_CMACRO1(x,y) \
        lua_pushnumber(L, y); \
        lua_setglobal(L, #x)

#define LUA_CTABLE1(x) \
        lua_newtable(L); \
        lua_setglobal(L, #x)

#define LUA_CADDCONST_VALUE_TO_TABLE1(const_name,const_value,namespace_name) \
        lua_getglobal(L,#namespace_name); \
        lua_pushstring(L,#const_name); \
        lua_pushnumber(L,const_value); \
        lua_settable(L,-3); \
        lua_pop(L,1)

void chi_mesh::lua_utils::RegisterLuaEntities(lua_State *L)
{
  //=================================== Unpartitioned Mesh
  LUA_FMACRO1(chiMeshCreateUnpartitioned1DOrthoMesh);
  LUA_FMACRO1(chiMeshCreateUnpartitioned2DOrthoMesh);
  LUA_FMACRO1(chiMeshCreateUnpartitioned3DOrthoMesh);

  unpartition_mesh_lua_utils::RegisterLuaEntities(L);

  //=================================== Logical Volume
  LUA_FMACRO1(chiLogicalVolumeCreate);
    LUA_CMACRO1(SPHERE       , 1);
    LUA_CMACRO1(SPHERE_ORIGIN, 2);
    LUA_CMACRO1(RPP          , 3);
    LUA_CMACRO1(RCC          , 4);
    LUA_CMACRO1(SURFACE      , 9);
    LUA_CMACRO1(BOOLEAN      , 10);
  LUA_FMACRO1(chiLogicalVolumePointSense);

  //=================================== Mesh handler
  LUA_FMACRO1(chiMeshHandlerCreate);
  LUA_FMACRO1(chiMeshHandlerSetCurrent);
  LUA_FMACRO1(chiMeshHandlerExportMeshToObj);
  LUA_FMACRO1(chiMeshHandlerExportMeshToVTK);

  //=================================== Surface Mesher
  LUA_FMACRO1(chiSurfaceMesherCreate);
    LUA_CMACRO1(SURFACEMESHER_PREDEFINED, 1);
  LUA_FMACRO1(chiSurfaceMesherExecute);
  LUA_FMACRO1(chiSurfaceMesherSetProperty);
    LUA_CMACRO1(MAX_AREA   , 1);

  //=================================== Volume Mesher
  LUA_FMACRO1(chiVolumeMesherCreate);
    LUA_CMACRO1(VOLUMEMESHER_EXTRUDER     , 4);
      LUA_CTABLE1(ExtruderTemplateType);
        LUA_CADDCONST_VALUE_TO_TABLE1(SURFACE_MESH,1,ExtruderTemplateType);
        LUA_CADDCONST_VALUE_TO_TABLE1(UNPARTITIONED_MESH,2,ExtruderTemplateType);
    LUA_CMACRO1(VOLUMEMESHER_UNPARTITIONED, 6);
  LUA_FMACRO1(chiVolumeMesherExecute);
  LUA_FMACRO1(chiVolumeMesherSetProperty);

    LUA_CMACRO1(FORCE_POLYGONS     , 1);
    LUA_CMACRO1(MESH_GLOBAL        , 2);
    LUA_CMACRO1(PARTITION_Z        , 3);
    LUA_CMACRO1(VOLUMEPARTITION_Y  , 4);
    LUA_CMACRO1(VOLUMEPARTITION_X  , 5);
    LUA_CMACRO1(CUTS_Z             , 6);
    LUA_CMACRO1(CUTS_Y             , 7);
    LUA_CMACRO1(CUTS_X             , 8);
    LUA_CMACRO1(PARTITION_TYPE     , 9);
      LUA_CMACRO1(KBA_STYLE_XYZ, 2);
      LUA_CMACRO1(PARMETIS     , 3);
    LUA_CMACRO1(EXTRUSION_LAYER    , 10);
    LUA_CMACRO1(MATID_FROMLOGICAL  , 11);
    LUA_CMACRO1(BNDRYID_FROMLOGICAL, 12);
    LUA_CMACRO1(MATID_FROM_LUA_FUNCTION  , 13);
    LUA_CMACRO1(BNDRYID_FROM_LUA_FUNCTION, 14);

  LUA_FMACRO1(chiVolumeMesherSetKBAPartitioningPxPyPz);
  LUA_FMACRO1(chiVolumeMesherSetKBACutsX);
  LUA_FMACRO1(chiVolumeMesherSetKBACutsY);
  LUA_FMACRO1(chiVolumeMesherSetKBACutsZ);

  LUA_FMACRO1(chiVolumeMesherSetMatIDToAll);
  LUA_FMACRO1(chiVolumeMesherSetupOrthogonalBoundaries);

    LUA_CTABLE1(OrthoBoundaryID);
      LUA_CADDCONST_VALUE_TO_TABLE1(XMAX,0,OrthoBoundaryID);
      LUA_CADDCONST_VALUE_TO_TABLE1(XMIN,1,OrthoBoundaryID);
      LUA_CADDCONST_VALUE_TO_TABLE1(YMAX,2,OrthoBoundaryID);
      LUA_CADDCONST_VALUE_TO_TABLE1(YMIN,3,OrthoBoundaryID);
      LUA_CADDCONST_VALUE_TO_TABLE1(ZMAX,4,OrthoBoundaryID);
      LUA_CADDCONST_VALUE_TO_TABLE1(ZMIN,5,OrthoBoundaryID);

  //=================================== Domain Decomposition
  LUA_FMACRO1(chiDecomposeSurfaceMeshPxPy);

  //=================================== Mesh Cutting
  LUA_FMACRO1(chiCutMesh);
  LUA_FMACRO1(chiCountMeshInLogicalVolume);

  //=================================== Field function interpolation
  ff_interpolation_lua_utils::RegisterLuaEntities(L);

  //=================================== Surface Mesh
  LUA_FMACRO1(chiSurfaceMeshCreate);
  LUA_FMACRO1(chiSurfaceMeshImportFromOBJFile);
  LUA_FMACRO1(chiSurfaceMeshImportFromTriangleFiles);

  LUA_FMACRO1(chiSurfaceMeshExportToObj);
  LUA_FMACRO1(chiSurfaceMeshExportPolyFile);

  LUA_FMACRO1(chiSurfaceMeshCheckCycles);
  LUA_FMACRO1(chiComputeLoadBalancing);


}