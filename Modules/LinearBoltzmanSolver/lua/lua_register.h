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
AddNamedConstantToNamespace(REFLECTING        ,3,LBSBoundaryTypes)
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
RegisterNamespace(LBSGroupset)
RegisterFunction(chiLBSCreateGroup)
RegisterFunction(chiLBSCreateGroupset)
RegisterFunction(chiLBSGroupsetAddGroups)
RegisterFunction(chiLBSGroupsetSetQuadrature)
RegisterFunction(chiLBSGroupsetSetAngleAggregationType)
AddNamedConstantToNamespace(ANGLE_AGG_SINGLE,1,LBSGroupset)
AddNamedConstantToNamespace(ANGLE_AGG_POLAR ,2,LBSGroupset)
RegisterFunction(chiLBSGroupsetSetAngleAggDiv)
RegisterFunction(chiLBSGroupsetSetGroupSubsets)
RegisterFunction(chiLBSGroupsetSetIterativeMethod)
RegisterFunction(chiLBSGroupsetSetResidualTolerance)
RegisterFunction(chiLBSGroupsetSetMaxIterations)
RegisterFunction(chiLBSGroupsetSetGMRESRestartIntvl)
RegisterFunction(chiLBSGroupsetSetEnableSweepLog)
RegisterFunction(chiLBSGroupsetSetWGDSA)
RegisterFunction(chiLBSGroupsetSetTGDSA)