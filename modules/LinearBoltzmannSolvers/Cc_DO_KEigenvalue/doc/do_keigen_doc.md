\defgroup LBSDOKeigen LBS Discrete Ordinates k-eigenvalue solver

# Discrete Ordinates k-eigenvalue solver

## 1 Available basic options

- "K_EIGEN_METHOD"
  - "power" Standard power iteration
  - "power1" Standard power iteration but with diffusion based synthetic acceleration,
             with the acceleration scheme using power iteration (PISA_PI)
  - "power2" Standard power iteration but with diffusion based synthetic acceleration,
             with the acceleration scheme using a non-linear k-eigenvalue formulation
             (PISA_NL)
  - "nonlinear" Pure non-linear k-eigenvalue formulation. Can be executed with
    Within-group DSA as well as Two-Grid acceleration.
  
\code
chiSolverSetBasicOption(phys1, "K_EIGEN_METHOD", "power")
chiSolverSetBasicOption(phys1, "K_EIGEN_METHOD", "power1")
chiSolverSetBasicOption(phys1, "K_EIGEN_METHOD", "power2")
chiSolverSetBasicOption(phys1, "K_EIGEN_METHOD", "nonlinear")
\endcode

- "K_EIGEN_RESET_SOLUTION". true/false \[Default:true\]. During a call to execute,
  this flag controls the initial setting for the flux moments solution. If set to
  true the flux moments will be initialized with a 1.0 for all the scalar moments. If
  set to false the solver will start by using the current `phi_old_local` in memory.
  This option is useful if you want to perform a few free power iterations before 
  using the non-linear method (which can be unstable without a proper initial guess).

  
### 1.1 Power iteration k-eigenvalue solver options
The following options are used when "K_EIGEN_METHOD" is set to "power", "power1" or
"power2".
- "PI_MAX_ITS". \[Default:100\]. Allowable number of power iterations.
- "PI_K_TOL". \[Default:1.0e-10\]. The relative change in `k_eff` below which to
  terminate.

\code
chiSolverSetBasicOption(phys1, "PI_MAX_ITS", 1000)
chiSolverSetBasicOption(phys1, "PI_K_TOL", 1.0e-10)
\endcode

#### Power Iteration Synthetic Acceleration (PISA) options
- "PISA_VERBOSE_LEVEL". \[Default:0\]. Level of verbosity for the acceleration 
  schemes. Generally, 0 does not print outer or inner iterations. 1 prints outer
  iterations only. 2 prints outer- and inner iterations.
- "PISA_MIP_L_ABS_TOL". \[Default:1.0e-10\] Tolerance for the iterative inversion
  of the diffusion operator.
- "PISA_MIP_L_MAX_ITS". \[Default:100\] Maximum iterations for the iterative inversion
  of the diffusion operator.

\code
chiSolverSetBasicOption(phys1, "PISA_VERBOSE_LEVEL", 2)
\endcode

**Using inner power iteration**
- "PISA_PI_K_TOL". \[Default:1.0e-10\] Iterative tolerance on the corrective 
  eigenvalue.
- "PISA_PI_MAX_ITS". \[Default:50\] Number of iterations to which the correction
  solve is limited to.

\code
chiSolverSetBasicOption(phys1, "PISA_PI_K_TOL", 1.0e-8)
chiSolverSetBasicOption(phys1, "PISA_PI_MAX_ITS", 10)
\endcode

**Using inner non-linear k-eigenvalue method**
- "PISA_NL_ABS_TOL", \[Default:1.0e-10\]. Non-linear absolute tolerance.
- "PISA_NL_REL_TOL", \[Default:1.0e-10\]. Non-linear relative tolerance.
- "PISA_NL_SOL_TOL", \[Default:1.0e-50\]. Non-linear solution tolerance.
- "PISA_NL_MAX_ITS", \[Default:50\]. Non-linear maximum iterations.
\n  

- "PISA_L_ABS_TOL", \[Default:1.0e-10\]. Non-linear inner linear absolute tolerance.
- "PISA_L_REL_TOL", \[Default:1.0e-10\]. Non-linear inner linear relative tolerance.
- "PISA_L_MAX_ITS", \[Default:50\]. Non-linear inner linear maximum iterations.


### 1.2 Non-Linear k-eigenvalue solver options
- "NLK_ABS_TOL" \[Default:1.0e-8\] Non-linear absolute tolerance.
- "NLK_REL_TOL" \[Default:1.0e-8\] Non-linear relative tolerance.
- "NLK_SOL_TOL" \[Default:1.0e-50\] Non-linear solution tolerance.
- "NLK_MAX_ITS" \[Default:50\] Non-linear maximum iterations.
\n  

- "NLK_L_REL_TOL" \[Default:1.0e-8\] Non-linear inner linear relative tolerance.
- "NLK_L_ABS_TOL" \[Default:1.0e-8\] Non-linear inner linear absolute tolerance.
- "NLK_L_DIV_TOL" \[Default:1.0e6\] Non-linear inner linear divergence tolerance.
- "NLK_L_MAX_ITS" \[Default:50\] Non-linear inner linear maximum iterations.
- "NLK_GMRES_RESTART_INTVL" \[Default:30\] Non-linear inner linear GMRes restart intvl.
- "NLK_GMRES_BRKDN_TOL" \[Default:1.0e6\] Non-linear inner linear GMRes breakdown tolerance.

\ingroup LuaLBS