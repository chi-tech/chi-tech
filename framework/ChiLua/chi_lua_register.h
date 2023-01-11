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
        lua_settable(this->consoleState,-3); \
        lua_pop(this->consoleState,1);

//#################################### MODULES
//module:MPI Utilities
//RegisterFunction(chiMPIBarrier)

//module:Logging Utilities
//RegisterFunction(chiLogSetVerbosity)
//RegisterFunction(chiLog)
//RegisterConstant(LOG_0,          1);
//RegisterConstant(LOG_0WARNING,   2);
//RegisterConstant(LOG_0ERROR,     3);
//RegisterConstant(LOG_0VERBOSE_0, 4);
//RegisterConstant(LOG_0VERBOSE_1, 5);
//RegisterConstant(ChiLog::LOG_LVL::LOG_0VERBOSE_2, 6);
//RegisterConstant(LOG_ALL,          7);
//RegisterConstant(LOG_ALLWARNING,   8);
//RegisterConstant(LOG_ALLERROR,     9);
//RegisterConstant(LOG_ALLVERBOSE_0, 10);
//RegisterConstant(LOG_ALLVERBOSE_1, 11);
//RegisterConstant(LOG_ALLVERBOSE_2, 12);

//module:Math Utilities
//RegisterFunction(chiLegendre)
//RegisterFunction(chiLegendreDerivative)
//RegisterFunction(chiYlm)

//module:Quadratures
//RegisterFunction(chiCreateLineQuadrature)
//RegisterFunction(chiCreateProductQuadrature)
//RegisterFunction(chiCreateCylindricalProductQuadrature)
//RegisterFunction(chiCreateSphericalProductQuadrature)
//    RegisterConstant(GAUSS_LEGENDRE,             1);
//    RegisterConstant(GAUSS_CHEBYSHEV,            2);
//    RegisterConstant(GAUSS_LEGENDRE_LEGENDRE,    3);
//    RegisterConstant(GAUSS_LEGENDRE_CHEBYSHEV,   4);
//    RegisterConstant(CUSTOM_QUADRATURE,          5);
//RegisterFunction(chiCreateCustomAngularQuadrature)
//RegisterFunction(chiGetProductQuadrature)
//RegisterFunction(chiCreateSLDFESQAngularQuadrature)
//RegisterFunction(chiLocallyRefineSLDFESQAngularQuadrature)
//RegisterFunction(chiPrintToPythonSLDFESQAngularQuadrature)
//RegisterFunction(chiOptimizeAngularQuadratureForPolarSymmetry)

//module:Orthogonal Meshes
//RegisterFunction(chiMeshCreateUnpartitioned1DOrthoMesh)
//RegisterFunction(chiMeshCreateUnpartitioned2DOrthoMesh)
//RegisterFunction(chiMeshCreateUnpartitioned3DOrthoMesh)

//module:Loading a Mesh from File
//RegisterFunction(chiUnpartitionedMeshFromVTU)
//RegisterFunction(chiUnpartitionedMeshFromEnsightGold)
//RegisterFunction(chiUnpartitionedMeshFromWavefrontOBJ)
//RegisterFunction(chiUnpartitionedMeshFromMshFormat)

//module:Manual Unpartitioned Mesh
//RegisterFunction(chiCreateEmptyUnpartitionedMesh)
//RegisterFunction(chiUnpartitionedMeshUploadVertex)
//RegisterFunction(chiUnpartitionedMeshUploadCell)
//RegisterFunction(chiUnpartitionedMeshFinalizeEmpty)
//RegisterFunction(chiDestroyUnpartitionedMesh)

//module:Logical Volumes
//    RegisterFunction(chiLogicalVolumeCreate)
//      RegisterConstant(SPHERE,   1);
//      RegisterConstant(SPHERE_ORIGIN,   2);
//      RegisterConstant(RPP,   3);
//      RegisterConstant(RCC,   4);
//      RegisterConstant(SURFACE,   9);
//      RegisterConstant(BOOLEAN,   10);

//module:Mesh Handler
//    RegisterFunction(chiMeshHandlerCreate)
//    RegisterFunction(chiMeshHandlerSetCurrent)
//    RegisterFunction(chiMeshHandlerExportMeshToObj)
//    RegisterFunction(chiMeshHandlerExportMeshToVTK)


//module:Surface Meshes
//    RegisterFunction(chiSurfaceMeshCreate)
//    RegisterFunction(chiSurfaceMeshImportFromOBJFile)
//    RegisterFunction(chiSurfaceMeshImportFromTriangleFiles)
//    RegisterFunction(chiSurfaceMeshExportToObj)
//    RegisterFunction(chiSurfaceMeshExportPolyFile)
//    RegisterFunction(chiSurfaceMeshCheckCycles)
//    RegisterFunction(chiComputeLoadBalancing)

//module:Meshers
//  SurfaceMesher
//    RegisterFunction(chiSurfaceMesherCreate)
//      RegisterConstant(SURFACEMESHER_PREDEFINED,   1);
//      RegisterConstant(SURFACEMESHER_DELAUNAY,   2);

//    RegisterFunction(chiSurfaceMesherExecute)
//    RegisterFunction(chiSurfaceMesherSetProperty)
//      RegisterConstant(MAX_AREA,   1);
//  VolumeMesher
//    RegisterFunction(chiVolumeMesherCreate)
//      RegisterConstant(VOLUMEMESHER_EXTRUDER,     4);
//      RegisterConstant(VOLUMEMESHER_UNPARTITIONED, 6)
//    RegisterNamespace(ExtruderTemplateType)
//      AddNamedConstantToNamespace(SURFACE_MESH,1,ExtruderTemplateType)
//      AddNamedConstantToNamespace(UNPARTITIONED_MESH,2,ExtruderTemplateType)
//
//    RegisterFunction(chiVolumeMesherExecute)
//    RegisterFunction(chiVolumeMesherSetProperty)
//      RegisterConstant(FORCE_POLYGONS,   1);
//      RegisterConstant(MESH_GLOBAL,   2);
//      RegisterConstant(PARTITION_Z,   3);
//      RegisterConstant(VOLUMEPARTITION_Y,   4);
//      RegisterConstant(VOLUMEPARTITION_X,   5);
//      RegisterConstant(CUTS_Z,   6);
//      RegisterConstant(CUTS_Y,   7);
//      RegisterConstant(CUTS_X,   8);
//      RegisterConstant(PARTITION_TYPE,   9);
//        RegisterConstant(KBA_STYLE_XYZ,   2);
//        RegisterConstant(PARMETIS,   3);
//      RegisterConstant(EXTRUSION_LAYER,   10);
//      RegisterConstant(MATID_FROMLOGICAL,   11);
//      RegisterConstant(BNDRYID_FROMLOGICAL, 12);
//    RegisterFunction(chiVolumeMesherSetKBAPartitioningPxPyPz)
//    RegisterFunction(chiVolumeMesherSetKBACutsX)
//    RegisterFunction(chiVolumeMesherSetKBACutsY)
//    RegisterFunction(chiVolumeMesherSetKBACutsZ)
//    RegisterFunction(chiVolumeMesherSetMatIDToAll)
//    RegisterFunction(chiVolumeMesherSetupOrthogonalBoundaries)
//    RegisterNamespace(OrthoBoundaryID)
//      AddNamedConstantToNamespace(XMAX,0,OrthoBoundaryID)
//      AddNamedConstantToNamespace(XMIN,1,OrthoBoundaryID)
//      AddNamedConstantToNamespace(YMAX,2,OrthoBoundaryID)
//      AddNamedConstantToNamespace(YMIN,3,OrthoBoundaryID)
//      AddNamedConstantToNamespace(ZMAX,4,OrthoBoundaryID)
//      AddNamedConstantToNamespace(ZMIN,5,OrthoBoundaryID)

//module:Experimental
//    RegisterFunction(chiDecomposeSurfaceMeshPxPy)
//    RegisterFunction(chiCutMesh)
//    RegisterFunction(chiCountMeshInLogicalVolume)

//module:Field-functions
//RegisterFunction(chiGetFieldFunctionHandleByName)
//RegisterFunction(chiGetFieldFunctionList)
//RegisterFunction(chiExportFieldFunctionToVTK)
//RegisterFunction(chiExportFieldFunctionToVTKG)
//RegisterFunction(chiExportMultiFieldFunctionToVTK)

//module:Field-function Manipulation
//    RegisterFunction(chiFFInterpolationCreate)
//      RegisterConstant(SLICE,   1);
//      RegisterConstant(LINE,   2);
//      RegisterConstant(VOLUME,   3);
//    RegisterFunction(chiFFInterpolationSetProperty)
//      RegisterConstant(ADD_FIELDFUNCTION,   0);
//      RegisterConstant(SLICE_POINT,   1);
//      RegisterConstant(SLICE_NORMAL,   2);
//      RegisterConstant(SLICE_TANGENT,   3);
//      RegisterConstant(SLICE_BINORM,   4);
//      RegisterConstant(OPERATION,   5);
//        RegisterConstant(OP_SUM,   10);
//        RegisterConstant(OP_AVG,   11);
//        RegisterConstant(OP_MAX,   12);
//        RegisterConstant(OP_SUM_LUA,   13);
//        RegisterConstant(OP_AVG_LUA,   14);
//        RegisterConstant(OP_MAX_LUA,   15);
//    RegisterConstant(LOGICAL_VOLUME,   8);
//    RegisterConstant(LINE_FIRSTPOINT,   11);
//    RegisterConstant(LINE_SECONDPOINT,   12);
//    RegisterConstant(LINE_NUMBEROFPOINTS,   13);
//    RegisterConstant(LINE_CUSTOM_ARRAY,   14);
//    RegisterFunction(chiFFInterpolationInitialize)
//    RegisterFunction(chiFFInterpolationExecute)
//    RegisterFunction(chiFFInterpolationExportPython)
//    RegisterFunction(chiFFInterpolationGetValue)

//module:Physics Utilities
//RegisterFunction(chiSolverInitialize)
//RegisterFunction(chiSolverExecute)
//RegisterFunction(chiSolverSetBasicOption)

//RegisterFunction(chiPhysicsAddMaterial)
//RegisterFunction(chiPhysicsMaterialAddProperty)
//RegisterFunction(chiPhysicsMaterialSetProperty)
//RegisterFunction(chiPhysicsMaterialGetProperty)
//RegisterFunction(chiPhysicsMaterialModifyTotalCrossSection)

//module:Transport interaction cross-sections
//\ref ChiXSFile Chi-Tech native cross-section file format
//RegisterFunction(chiPhysicsTransportXSCreate)
//RegisterFunction(chiPhysicsTransportXSSet)
//RegisterFunction(chiPhysicsTransportXSMakeCombined)
//RegisterFunction(chiPhysicsTransportXSGet)
//RegisterFunction(chiPhysicsTransportXSExportToChiTechFormat)
//
////Property indices
//RegisterConstant(SCALAR_VALUE,           1);
//RegisterConstant(TRANSPORT_XSECTIONS,    10);
//RegisterConstant(ISOTROPIC_MG_SOURCE,    11);
//
////Operation indices
//RegisterConstant(SINGLE_VALUE,            0);
//RegisterConstant(FROM_ARRAY,              1);
//RegisterConstant(SIMPLEXS0,              20);
//RegisterConstant(SIMPLEXS1,              21);
//RegisterConstant(EXISTING,               22);
//RegisterConstant(CHI_XSFILE,             23);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "../../modules/module_lua_register.h"

#endif
//module:Test scripts
//RegisterFunction(chiLuaTest)
//
//
//RegisterNamespace(LuaNamespace)
//AddNamedConstantToNamespace(Name,1,LuaNamespace)

//string: This file was generated using the script: CHI_TECH/ChiLua/chi_lua_docbuild.lua
