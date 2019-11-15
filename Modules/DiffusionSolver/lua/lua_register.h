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