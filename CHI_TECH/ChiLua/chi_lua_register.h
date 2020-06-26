#define RegisterFunction(x) \
        int x(lua_State *L); \
        lua_register(this->consoleState, #x, x);

#define RegisterConstant(x,y) \
        lua_pushnumber(this->consoleState,y); \
        lua_setglobal(this->consoleState, #x);

#define RegisterNamespace(x) \
        lua_newtable(this->consoleState); \
        lua_setglobal(this->consoleState, #x);

#define AddNamedConstantToNamespace(const_name,const_value,namespace_name) \
        lua_getglobal(this->consoleState,#namespace_name); \
        lua_pushstring(this->consoleState,#const_name); \
        lua_pushnumber(this->consoleState,const_value); \
        lua_settable(this->consoleState,-3);




//#################################### MODULES
//module:Math Utilities
RegisterFunction(chiLegendre)
RegisterFunction(chiLegendreDerivative)
RegisterFunction(chiYlm)
RegisterFunction(chiCreateQuadrature)
RegisterFunction(chiCreateProductQuadrature)
    RegisterConstant(GAUSS_LEGENDRE,             1);
    RegisterConstant(GAUSS_CHEBYSHEV,            2);
    RegisterConstant(GAUSS_LEGENDRE_LEGENDRE,    3);
    RegisterConstant(GAUSS_LEGENDRE_CHEBYSHEV,   4);
    RegisterConstant(CUSTOM_QUADRATURE,          5);
RegisterFunction(chiCreateCustomAngularQuadrature)
RegisterFunction(chiGetProductQuadrature)

//module:Mesh Macros
RegisterFunction(chiMeshCreate1DSlabMesh)
RegisterFunction(chiMeshCreate2DOrthoMesh)
RegisterFunction(chiMeshCreate3DOrthoMesh)
RegisterFunction(chiUnpartitionedMeshFromVTU)
RegisterFunction(chiUnpartitionedMeshFromEnsightGold)


//module:Mesh Utilities
RegisterFunction(chiEdgeLoopSplitByAngle)
//  LineMesh
    RegisterFunction(chiLineMeshCreateFromLoop)
    RegisterFunction(chiLineMeshCreateFromArray)
//  Logical volumes
    RegisterFunction(chiLogicalVolumeCreate)
      RegisterConstant(SPHERE,   1);
      RegisterConstant(SPHERE_ORIGIN,   2);
      RegisterConstant(RPP,   3);
      RegisterConstant(RCC,   4);
      RegisterConstant(SURFACE,   9);
      RegisterConstant(BOOLEAN,   10);
//  Handler
    RegisterFunction(chiMeshHandlerCreate)
    RegisterFunction(chiMeshHandlerGetSurfaceFromCollection)
    RegisterFunction(chiMeshHandlerSetCurrent)
//  Region
    RegisterFunction(chiRegionCreate)
    RegisterFunction(chiRegionAddSurfaceBoundary)
    RegisterFunction(chiRegionAddLineBoundary)
    RegisterFunction(chiRegionAddEmptyBoundary)
    RegisterFunction(chiRegionGetBoundarySurfaceMesh)
    RegisterFunction(chiRegionExportMeshToPython)
    RegisterFunction(chiRegionExportMeshToObj)
    RegisterFunction(chiRegionExportMeshToVTK)
//  SurfaceMesh
    RegisterFunction(chiSurfaceMeshCreate)
    RegisterFunction(chiSurfaceMeshCreateFromArrays)
    RegisterFunction(chiSurfaceMeshImportFromOBJFile)
    RegisterFunction(chiSurfaceMeshImportFromTriangleFiles)
    RegisterFunction(chiSurfaceMeshImportFromMshFiles)
    RegisterFunction(chiSurfaceMeshExportToObj)
    RegisterFunction(chiSurfaceMeshExportPolyFile)
    RegisterFunction(chiSurfaceMeshGetEdgeLoops)
    RegisterFunction(chiSurfaceMeshGetEdgeLoopsPoly)
    RegisterFunction(chiSurfaceMeshSplitByPatch)
    RegisterFunction(chiSurfaceMeshExtractOpenEdgesToObj)
    RegisterFunction(chiSurfaceMeshCheckCycles)
    RegisterFunction(chiComputeLoadBalancing)
//  SurfaceMesher
    RegisterFunction(chiSurfaceMesherCreate)
      RegisterConstant(SURFACEMESHER_PREDEFINED,   1);
      RegisterConstant(SURFACEMESHER_DELAUNAY,   2);
      RegisterConstant(SURFACEMESHER_TRIANGLE,   3);
    RegisterFunction(chiSurfaceMesherExecute)
    RegisterFunction(chiSurfaceMesherSetProperty)
      RegisterConstant(MAX_AREA,   1);
      RegisterConstant(PARTITION_X,   2);
      RegisterConstant(PARTITION_Y,   3);
      RegisterConstant(CUT_X,   4);
      RegisterConstant(CUT_Y,   5);
    RegisterFunction(chiSurfaceMesherExportToObj)
//  VolumeMesher
    RegisterFunction(chiVolumeMesherCreate)
      RegisterConstant(VOLUMEMESHER_LINEMESH1D,   1);
      RegisterConstant(VOLUMEMESHER_PREDEFINED2D, 3);
      RegisterConstant(VOLUMEMESHER_EXTRUDER,     4);
      RegisterConstant(VOLUMEMESHER_PREDEFINED3D, 5)
    RegisterFunction(chiVolumeMesherExecute)
    RegisterFunction(chiVolumeMesherSetProperty)
      RegisterConstant(FORCE_POLYGONS,   1);
      RegisterConstant(MESH_GLOBAL,   2);
      RegisterConstant(PARTITION_Z,   3);
      RegisterConstant(VOLUMEPARTITION_Y,   4);
      RegisterConstant(VOLUMEPARTITION_X,   5);
      RegisterConstant(CUTS_Z,   6);
      RegisterConstant(CUTS_Y,   7);
      RegisterConstant(CUTS_X,   8);
      RegisterConstant(PARTITION_TYPE,   9);
        RegisterConstant(KBA_STYLE_XY,   1);
        RegisterConstant(KBA_STYLE_XYZ,   2);
        RegisterConstant(PARMETIS,   3);
      RegisterConstant(EXTRUSION_LAYER,   10);
      RegisterConstant(MATID_FROMLOGICAL,   11);
      RegisterConstant(BNDRYID_FROMLOGICAL, 12);
//  Domain Decomposition
    RegisterFunction(chiDomDecompose2D)
    RegisterFunction(chiDecomposeSurfaceMeshPxPy)
//module:Field-function Manipulation
    RegisterFunction(chiFFInterpolationCreate)
      RegisterConstant(SLICE,   1);
      RegisterConstant(LINE,   2);
      RegisterConstant(VOLUME,   3);
    RegisterFunction(chiFFInterpolationSetProperty)
      RegisterConstant(ADD_FIELDFUNCTION,   0);
      RegisterConstant(SLICE_POINT,   1);
      RegisterConstant(SLICE_NORMAL,   2);
      RegisterConstant(SLICE_TANGENT,   3);
      RegisterConstant(SLICE_BINORM,   4);
      RegisterConstant(OPERATION,   5);
        RegisterConstant(OP_SUM,   10);
        RegisterConstant(OP_AVG,   11);
        RegisterConstant(OP_MAX,   12);
        RegisterConstant(OP_SUM_LUA,   13);
        RegisterConstant(OP_AVG_LUA,   14);
        RegisterConstant(OP_MAX_LUA,   15);
    RegisterConstant(LOGICAL_VOLUME,   8);
    RegisterConstant(LINE_FIRSTPOINT,   11);
    RegisterConstant(LINE_SECONDPOINT,   12);
    RegisterConstant(LINE_NUMBEROFPOINTS,   13);
    RegisterConstant(LINE_CUSTOM_ARRAY,   14);
    RegisterFunction(chiFFInterpolationInitialize)
    RegisterFunction(chiFFInterpolationExecute)
    RegisterFunction(chiFFInterpolationExportPython)
    RegisterFunction(chiFFInterpolationGetValue)


//module:MPI Utilities
RegisterFunction(chiMPIBarrier)

//module:Logging Utilities
RegisterFunction(chiLogSetVerbosity)
RegisterFunction(chiLog)
RegisterConstant(LOG_0,          1);
RegisterConstant(LOG_0WARNING,   2);
RegisterConstant(LOG_0ERROR,     3);
RegisterConstant(LOG_0VERBOSE_0, 4);
RegisterConstant(LOG_0VERBOSE_1, 5);
RegisterConstant(LOG_0VERBOSE_2, 6);
RegisterConstant(LOG_ALL,          7);
RegisterConstant(LOG_ALLWARNING,   8);
RegisterConstant(LOG_ALLERROR,     9);
RegisterConstant(LOG_ALLVERBOSE_0, 10);
RegisterConstant(LOG_ALLVERBOSE_1, 11);
RegisterConstant(LOG_ALLVERBOSE_2, 12);

//module:Physics Utilities
RegisterFunction(chiSolverAddRegion)
RegisterFunction(chiSolverExecute)
RegisterFunction(chiPhysicsAddMaterial)
RegisterFunction(chiPhysicsMaterialAddProperty)
RegisterFunction(chiPhysicsMaterialSetProperty)
RegisterFunction(chiPhysicsMaterialGetProperty)
RegisterFunction(chiPhysicsTransportXSCreate)
RegisterFunction(chiPhysicsTransportXSSet)
RegisterFunction(chiPhysicsTransportXSMakeCombined)
RegisterFunction(chiGetFieldFunctionList)
RegisterFunction(chiExportFieldFunctionToVTK)
RegisterFunction(chiExportFieldFunctionToVTKG)
RegisterFunction(chiExportMultiFieldFunctionToVTKG)


//Property indices
RegisterConstant(SCALAR_VALUE,           1);
RegisterConstant(TRANSPORT_XSECTIONS,    10);
RegisterConstant(ISOTROPIC_MG_SOURCE,    11);

//Operation indices
RegisterConstant(SINGLE_VALUE,            0);
RegisterConstant(FROM_ARRAY,              1);
RegisterConstant(SIMPLEXS0,              20);
RegisterConstant(SIMPLEXS1,              21);
RegisterConstant(PDT_XSFILE,             22);
RegisterConstant(EXISTING,               23);
RegisterConstant(CHI_XSFILE,             24);

//Unit tests
RegisterFunction(chiExecuteUnitTests)

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "../../Modules/module_lua_register.h"

#endif
//module:Test scripts
RegisterFunction(chiLuaTest)


RegisterNamespace(LuaNamespace)
AddNamedConstantToNamespace(Name,1,LuaNamespace)

