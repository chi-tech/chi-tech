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
    RegisterConstant(GAUSS_LEGENDRE,   1);
    RegisterConstant(GAUSS_CHEBYSHEV,   2);
RegisterFunction(chiCreateProductQuadrature)
    RegisterConstant(GAUSS_LEGENDRE,             1);
    RegisterConstant(GAUSS_CHEBYSHEV,            2);
    RegisterConstant(GAUSS_LEGENDRE_LEGENDRE,    3);
    RegisterConstant(GAUSS_LEGENDRE_CHEBYSHEV,   4);


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
    RegisterFunction(chiRegionGetBoundarySurfaceMesh)
    RegisterFunction(chiRegionExportMeshToPython)
    RegisterFunction(chiRegionExportMeshToObj)
    RegisterFunction(chiRegionExportMeshToVTK)
//  SurfaceMesh
    RegisterFunction(chiSurfaceMeshCreate)
    RegisterFunction(chiSurfaceMeshImportFromOBJFile)
    RegisterFunction(chiSurfaceMeshImportFromTriangleFiles)
    RegisterFunction(chiSurfaceMeshExportToObj)
    RegisterFunction(chiSurfaceMeshExportPolyFile)
    RegisterFunction(chiSurfaceMeshGetEdgeLoops)
    RegisterFunction(chiSurfaceMeshGetEdgeLoopsPoly)
    RegisterFunction(chiSurfaceMeshSplitByPatch)
    RegisterFunction(chiSurfaceMeshExtractOpenEdgesToObj)
    RegisterFunction(chiSurfaceMeshCheckCycles)
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
      RegisterConstant(VOLUMEMESHER_PREDEFINED2D,   3);
      RegisterConstant(VOLUMEMESHER_EXTRUDER,   4);
    RegisterFunction(chiVolumeMesherExecute)
    RegisterFunction(chiVolumeMesherSetProperty)
      RegisterConstant(FORCE_POLYGONS,   1);
      RegisterConstant(MESH_GLOBAL,   2);
      RegisterConstant(PARTITION_Z,   3);
      RegisterConstant(EXTRUSION_LAYER,   10);
      RegisterConstant(MATID_FROMLOGICAL,   11);
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
    RegisterConstant(LOGICAL_VOLUME,   8);
    RegisterConstant(LINE_FIRSTPOINT,   11);
    RegisterConstant(LINE_SECONDPOINT,   12);
    RegisterConstant(LINE_NUMBEROFPOINTS,   13);
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
RegisterFunction(chiSolverAddFieldFunction)
RegisterFunction(chiGetFieldFunctionList)
RegisterFunction(chiExportFieldFunctionToVTK)
RegisterFunction(chiExportFieldFunctionToVTKG)
RegisterFunction(chiExportMultiFieldFunctionToVTKG)


//Property indices
RegisterConstant(THERMAL_CONDUCTIVITY,   0);
RegisterConstant(SCALAR_VALUE,           1);
RegisterConstant(TRANSPORT_XSECTIONS,    10);
RegisterConstant(ISOTROPIC_MG_SOURCE,    11);

//Operation indices
RegisterConstant(SINGLE_VALUE,            0);
RegisterConstant(FROM_ARRAY,              1);
RegisterConstant(SIMPLEXS0,              20);
RegisterConstant(SIMPLEXS1,              21);
RegisterConstant(PDT_XSFILE,             22);


//module:Monte-carlo Utilities
RegisterFunction(chiMonteCarlonCreateSolver)
RegisterFunction(chiMonteCarlonCreateSource)
RegisterFunction(chiMonteCarlonInitialize)
RegisterFunction(chiMonteCarlonExecute)
  RegisterConstant(MC_POINT_SRC,             1);
  RegisterConstant(MC_BNDRY_SRC,             2);
    RegisterConstant(MC_ALL_BOUNDARIES,         -1);
  RegisterConstant(MC_LOGICVOL_SRC,          3);
  RegisterConstant(MC_RESID_SRC,             4);
RegisterFunction(chiMonteCarlonSetProperty)
  RegisterConstant(MC_NUM_PARTICLES,             1);
  RegisterConstant(MC_TFC_UPDATE_INTVL,          2);
  RegisterConstant(MC_MONOENERGETIC,             3);
  RegisterConstant(MC_SCATTERING_ORDER,          4);
  RegisterConstant(MC_FORCE_ISOTROPIC,           5);
  RegisterConstant(MC_GROUP_BOUNDS,              6);
  RegisterConstant(MC_TALLY_MERGE_INTVL,         7);
  RegisterConstant(MC_TALLY_MULTIPLICATION_FACTOR, 8);

//module:Diffusion
RegisterFunction(chiDiffusionCreateSolver)
RegisterFunction(chiDiffusionInitialize)
RegisterFunction(chiDiffusionExecute)
RegisterFunction(chiDiffusionSetProperty)
  RegisterConstant(DISCRETIZATION_METHOD,   1);
    RegisterConstant(PWLC,      3);
    RegisterConstant(PWLD_MIP,   4);
  RegisterConstant(MAX_ITERS,               2);
  RegisterConstant(RESIDUAL_TOL,            3);
  RegisterConstant(BOUNDARY_TYPE,           4);
    RegisterConstant(DIFFUSION_REFLECTING,-1);
    RegisterConstant(DIFFUSION_DIRICHLET, -2);
    RegisterConstant(DIFFUSION_NEUMANN,   -3);
    RegisterConstant(DIFFUSION_VACUUM,    -4);
    RegisterConstant(DIFFUSION_ROBIN,     -5);
  RegisterConstant(PROPERTY_D_MAP,          5);
  RegisterConstant(PROPERTY_Q_MAP,          6);
  RegisterConstant(PROPERTY_SIGMAA_MAP,     7);

//module:Linear Boltzman Solver
RegisterFunction(chiLBSCreateSolver)

RegisterFunction(chiLBSSetProperty)
RegisterConstant(DISCRETIZATION_METHOD,   1);
  RegisterConstant(PWLD1D,   4);
  RegisterConstant(PWLD2D,   5);
  RegisterConstant(PWLD3D,   6);
RegisterConstant(PARTITION_METHOD,   2);
  RegisterConstant(SERIAL,   1);
  RegisterConstant(FROM_SURFACE,   2);
RegisterConstant(BOUNDARY_CONDITION,   3);
  RegisterConstant(XMAX,   31);
  RegisterConstant(XMIN,   32);
  RegisterConstant(YMAX,   33);
  RegisterConstant(YMIN,   34);
  RegisterConstant(ZMAX,   35);
  RegisterConstant(ZMIN,   36);


  RegisterNamespace(LBSBoundaryTypes);
  AddNamedConstantToNamespace(VACUUM            ,1,LBSBoundaryTypes)
  AddNamedConstantToNamespace(INCIDENT_ISOTROPIC,2,LBSBoundaryTypes)
//
//    RegisterConstant(VACUUM,               301);
//    RegisterConstant(INCIDENT_ISOTROPIC,   302);

RegisterConstant(GROUPSET_ITERATIVEMETHOD,   101);
  RegisterConstant(NPT_CLASSICRICHARDSON,          1);
  RegisterConstant(NPT_CLASSICRICHARDSON_CYCLES,   2);
  RegisterConstant(NPT_GMRES,                      3);
  RegisterConstant(NPT_GMRES_CYCLES,               4);
RegisterConstant(GROUPSET_TOLERANCE,   102);
RegisterConstant(GROUPSET_MAXITERATIONS,   103);
RegisterConstant(GROUPSET_GMRESRESTART_INTVL,   104);
RegisterConstant(GROUPSET_SUBSETS,105);
RegisterConstant(GROUPSET_WGDSA,106);
RegisterConstant(GROUPSET_TGDSA,107);
RegisterConstant(GROUPSET_WGDSA_MAXITERATIONS,108);
RegisterConstant(GROUPSET_TGDSA_MAXITERATIONS,109);
RegisterConstant(GROUPSET_WGDSA_TOLERANCE,110);
RegisterConstant(GROUPSET_TGDSA_TOLERANCE,111);

RegisterConstant(SCATTERING_ORDER,   4);
RegisterConstant(SWEEP_EAGER_LIMIT,   5);
RegisterFunction(chiLBSInitialize)
RegisterFunction(chiLBSExecute)
RegisterFunction(chiLBSGetFieldFunctionList)
RegisterFunction(chiLBSGetScalarFieldFunctionList)

//module:Linear Boltzman Solver - Groupset manipulation
//\ref LuaLBSGroupsets Main page
RegisterFunction(chiLBSCreateGroup)
RegisterFunction(chiLBSCreateGroupset)
RegisterFunction(chiLBSGroupsetAddGroups)
RegisterFunction(chiLBSGroupsetSetQuadrature)
RegisterFunction(chiLBSGroupsetSetAngleAggDiv)
RegisterFunction(chiLBSGroupsetSetGroupSubsets)
RegisterFunction(chiLBSGroupsetSetIterativeMethod)
RegisterFunction(chiLBSGroupsetSetResidualTolerance)
RegisterFunction(chiLBSGroupsetSetMaxIterations)
RegisterFunction(chiLBSGroupsetSetGMRESRestartIntvl)
RegisterFunction(chiLBSGroupsetSetWGDSA)
RegisterFunction(chiLBSGroupsetSetTGDSA)

//module:Test scripts
RegisterFunction(chiLuaTest)


RegisterNamespace(LuaNamespace)
AddNamedConstantToNamespace(Name,1,LuaNamespace)

