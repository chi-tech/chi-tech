#include "../../../CHI_LUA/chi_lua.h"

#include "CHI_MODULES/LinearBoltzmanSolver/lbs_linear_boltzman_solver.h"
#include "../../../CHI_PHYSICS/chi_physics.h"
#include "../../../ChiMath/chi_math.h"
#include <chi_log.h>

extern CHI_PHYSICS chi_physics_handler;
extern ChiMath    chi_math_handler;
extern CHI_LOG     chi_log;

/** \defgroup LuaLBSGroupsets LBS Groupsets

The code below is an example of a complete specification of a groupset.

\code
--===================================== Setup physics
phys1 = chiLBSransportCreateSolver()
chiSolverAddRegion(phys1,region1)

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
chiLBSSetProperty(phys1,SCATTERING_ORDER,1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)
pquad1 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,8, 8)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)

cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,15)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad0)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-4)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)
chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")
\endcode

Groupsets segregate the code into pieces arranged by the number of groups
it contains. A great deal of care must be taken with intergroupset transfer
since the order in which the groupsets are executed determine what information
will be available to them.

\ingroup LuaNPT*/

//###################################################################
/**Create a groupset.
\param SolverIndex int Handle to the solver for which the set is to be created.

##_

Example:
\code
gs0 = chiLBSCreateGroupset(phys1)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiLBSCreateGroupset(lua_State *L)
{
  int solver_index = lua_tonumber(L,1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      fprintf(stderr,"ERROR: Incorrect solver-type"
                     "in chiLBSCreateGroupset\n");
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiLBSCreateGroupset\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Create groupset
  LBS_GROUPSET* newgs = new LBS_GROUPSET;
  solver->group_sets.push_back(newgs);

  lua_pushnumber(L,solver->group_sets.size()-1);
  return 1;
}



//###################################################################
/**Create a group.
\param SolverIndex int Handle to the solver for which the group
is to be created.

##_

Example:
\code
grp[g] = chiLBSCreateGroup(phys1)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiLBSCreateGroup(lua_State *L)
{
  int solver_index = lua_tonumber(L,1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      fprintf(stderr,"ERROR: Incorrect solver-type"
                     "in chiLBSCreateGroup\n");
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiLBSCreateGroup\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Create groupset
  LBS_GROUP* newgs = new LBS_GROUP;
  solver->groups.push_back(newgs);
  newgs->id = solver->groups.size()-1;
  lua_pushnumber(L,newgs->id);
  return 1;
}





//###################################################################
/**Adds a block of groups to a groupset.
\param SolverIndex int Handle to the solver for which the group
is to be created.
\param GroupsetIndex int Handle to the groupset to which the group is
 to be added.
\param FromIndex int From which group.
\param ToIndex int To which group.

##_

Example:
\code
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

chiLBSGroupsetAddGroups(phys1,cur_gs,0,15)
\endcode


\ingroup LuaLBSGroupsets
*/
int chiLBSGroupsetAddGroups(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 4)
  {
    fprintf(stderr,"ERROR: Invalid amount of arguments, %d,"
                   "in chiLBSGroupsetAddGroups"
                   " (4 required)\n",num_args);
    exit(EXIT_FAILURE);
  }
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int from = lua_tonumber(L,3);
  int to   = lua_tonumber(L,4);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      fprintf(stderr,"ERROR: Incorrect solver-type in "
                     "chiLBSGroupsetAddGroups\n");
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiLBSGroupsetAddGroups\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  LBS_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to groupset"
                   "in chiLBSGroupsetAddGroups\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Add the groups
  if (to<from)
  {
    chi_log.Log(LOG_0ERROR)
    << "No groups added to groupset in chiLBSGroupsetAddGroups. "
       "This is triggered when groups are added with the \"to\" "
       "field being less than the \"from\" field.";
    exit(EXIT_FAILURE);
  }


  for (unsigned k=from; k<=to; k++)
  {
    LBS_GROUP* group;
    //================================= Check valid group
    try {
      group = solver->groups.at(k);
    }
    catch (std::out_of_range o)
    {
      fprintf(stderr,"ERROR: Invalid group added to groupset"
                     "in chiLBSGroupsetAddGroups\n");
      exit(EXIT_FAILURE);
    }

    groupset->groups.push_back(group);
  }
  return 0;
}




//###################################################################
/**Sets the product quadrature used for the groupset
\param SolverIndex int Handle to the solver for which the group
is to be created.
\param GroupsetIndex int Handle to the groupset to which the group is
 to be added.
\param QuadratureIndex int Handle to the quadrature to be set for this
 groupset.


##_

Example:
\code
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)

chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad0)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiLBSGroupsetSetQuadrature(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiLBSGroupsetSetQuadrature",3,num_args);

  LuaCheckNilValue("chiLBSGroupsetSetQuadrature",L,1);
  LuaCheckNilValue("chiLBSGroupsetSetQuadrature",L,2);
  LuaCheckNilValue("chiLBSGroupsetSetQuadrature",L,3);


  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int prquad_index = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      fprintf(stderr,"ERROR: Incorrect solver-type"
                     "in chiLBSGroupsetSetQuadrature\n");
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiLBSGroupsetSetQuadrature\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  LBS_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to groupset"
                   "in chiLBSGroupsetSetQuadrature\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to quadrature
  chi_math::ProductQuadrature* prodquad;
  try{
    prodquad = chi_math_handler.product_quadratures.at(prquad_index);
  }
  catch (std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to Product Quadrature"
                   "in chiLBSGroupsetSetQuadrature\n");
    exit(EXIT_FAILURE);
  }

  groupset->quadrature = prodquad;

  chi_log.Log(LOG_0)
    << "Groupset " << grpset_index
    << " quadrature set to quadrature with "
    << prodquad->azimu_ang.size()
    << " azimuthal angles and "
    << prodquad->polar_ang.size()
    << " polar angles. ";


  return 0;
}

//###################################################################
/**Sets the angle aggregation divisions
\param SolverIndex int Handle to the solver for which the group
is to be created.
\param GroupsetIndex int Handle to the groupset to which the group is
 to be added.
\param NumDiv int Number of divisions to use for the angle aggregation.

Note: by default polar aggregation will combine all polar angles in a hemisphere
 for a given azimuthal angleset. Therefore if there are 24 polar angles and
 4 azimuthal angles the default polar aggregation will create 8 anglesets
 (2 per quadrant to allow top and bottom hemisphere) and each angleset will have the
 12 polar angles associated with a hemisphere. When the number of divisions is
 greater than 1 then the polar angles will be split into divisions. For example
 if the number of divisions is 2 then more angleset will be created, this time
 having 6 polar angles per angleset.

##_

Example:
\code
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiLBSGroupsetSetAngleAggDiv(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiLBSGroupsetSetAngleAggDiv",3,num_args);

  LuaCheckNilValue("chiLBSGroupsetSetAngleAggDiv",L,1);
  LuaCheckNilValue("chiLBSGroupsetSetAngleAggDiv",L,2);
  LuaCheckNilValue("chiLBSGroupsetSetAngleAggDiv",L,3);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int num_div = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiLBSGroupsetSetAngleAggDiv";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiLBSGroupsetSetAngleAggDiv";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  LBS_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiLBSGroupsetSetAngleAggDiv";
    exit(EXIT_FAILURE);
  }

  //============================================= Bounds checking
  if (num_div <= 0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid number of divisions "
      << "in call to chiLBSGroupsetSetAngleAggDiv. Must be >= 1.";
  }

  groupset->master_num_ang_subsets = num_div;

  chi_log.Log(LOG_0)
    << "Groupset " << grpset_index << " angle aggregation divisions "
    << "set to " << num_div;

  return 0;
}

//###################################################################
/**Sets the number of group-subsets to use for groupset. Default 1.
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param NumDiv int Number of divisions of the groupset to use.

##_

Example:
\code
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiLBSGroupsetSetGroupSubsets(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiLBSGroupsetSetGroupSubsets",3,num_args);

  LuaCheckNilValue("chiLBSGroupsetSetGroupSubsets",L,1);
  LuaCheckNilValue("chiLBSGroupsetSetGroupSubsets",L,2);
  LuaCheckNilValue("chiLBSGroupsetSetGroupSubsets",L,3);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int num_div = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiLBSGroupsetSetGroupSubsets";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiLBSGroupsetSetGroupSubsets";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  LBS_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiLBSGroupsetSetGroupSubsets";
    exit(EXIT_FAILURE);
  }

  //============================================= Bounds checking
  if (num_div <= 0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid number of subsets "
      << "in call to chiLBSGroupsetSetGroupSubsets. Must be >= 1.";
  }

  groupset->master_num_grp_subsets = num_div;

  chi_log.Log(LOG_0)
    << "Groupset " << grpset_index << " subset divisions "
    << "set to " << num_div;

  return 0;
}

//###################################################################
/**Sets the number of group-subsets to use for groupset. Default 1.
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param IterativeMethod int Iteritve method identifier.

##_

### IterativeMethod
NPT_CLASSICRICHARDSON\n
Standard source iteration.\n\n

NPT_GMRES\n
Generalized Minimal Residual formulation for iterations.\n\n

Example:
\code
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_CLASSICRICHARDSON)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiLBSGroupsetSetIterativeMethod(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiLBSGroupsetSetIterativeMethod",3,num_args);

  LuaCheckNilValue("chiLBSGroupsetSetIterativeMethod",L,1);
  LuaCheckNilValue("chiLBSGroupsetSetIterativeMethod",L,2);
  LuaCheckNilValue("chiLBSGroupsetSetIterativeMethod",L,3);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int iter_method  = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiLBSGroupsetSetGroupSubsets";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiLBSGroupsetSetGroupSubsets";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  LBS_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiLBSGroupsetSetGroupSubsets";
    exit(EXIT_FAILURE);
  }

  if (iter_method == NPT_CLASSICRICHARDSON)
  {
    groupset->iterative_method = NPT_CLASSICRICHARDSON;
  }
  else if (iter_method == NPT_GMRES)
  {
    groupset->iterative_method = NPT_GMRES;
  }
  else
  {
    chi_log.Log(LOG_0ERROR)
      << "Unsupported iterative method specified in call to "
      << "chiLBSGroupsetSetIterativeMethod.";
    exit(EXIT_FAILURE);
  }

  return 0;
}

//###################################################################
/**Sets the residual tolerance for the iterative method of the groupset.
 *
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param ResidualTol float residual tolerance (default 1.0e-6)

Note this tolerance also gets used for classic-richardson pointwise convergence
tolerance.

##_

Example:
\code
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-4)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiLBSGroupsetSetResidualTolerance(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiLBSGroupsetSetResidualTolerance",3,num_args);

  LuaCheckNilValue("chiLBSGroupsetSetResidualTolerance",L,1);
  LuaCheckNilValue("chiLBSGroupsetSetResidualTolerance",L,2);
  LuaCheckNilValue("chiLBSGroupsetSetResidualTolerance",L,3);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  double resid_tol = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiLBSGroupsetSetGroupSubsets";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiLBSGroupsetSetGroupSubsets";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  LBS_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiLBSGroupsetSetGroupSubsets";
    exit(EXIT_FAILURE);
  }

  //============================================= Bounds checking
  if (resid_tol < 0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid residual tolerance specified. Must be greater >= 0.0";
    exit(EXIT_FAILURE);
  }

  groupset->residual_tolerance = resid_tol;

  char buff[100];
  sprintf(buff,"%.4e",resid_tol);

  chi_log.Log(LOG_0)
    << "Groupset " << grpset_index << " residual tolerance "
    << "set to " << buff;

  return 0;
}

//###################################################################
/**Sets the maximum number of iterations for the groupset iterative method.
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param Numiter int Maximum number of iterations. Default 1000.

##_

Example:
\code
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,200)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiLBSGroupsetSetMaxIterations(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiLBSGroupsetSetMaxIterations",3,num_args);

  LuaCheckNilValue("chiLBSGroupsetSetMaxIterations",L,1);
  LuaCheckNilValue("chiLBSGroupsetSetMaxIterations",L,2);
  LuaCheckNilValue("chiLBSGroupsetSetMaxIterations",L,3);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int num_iter = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiLBSGroupsetSetMaxIterations";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiLBSGroupsetSetMaxIterations";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  LBS_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiLBSGroupsetSetMaxIterations";
    exit(EXIT_FAILURE);
  }

  //============================================= Bounds checking
  if (num_iter < 0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid number of iterations "
      << "in call to chiLBSGroupsetSetMaxIterations. Must be >= 0.";
  }

  groupset->max_iterations = num_iter;

  chi_log.Log(LOG_0)
    << "Groupset " << grpset_index << " max # iterations "
    << "set to " << num_iter;

  return 0;
}

//###################################################################
/**Sets the restart interval for GMRES if applied to the groupset.
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param Intvl int Interval to use for GMRES restarts. Default 30.

##_

Example:
\code
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,15)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiLBSGroupsetSetGMRESRestartIntvl(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiLBSGroupsetSetGMRESRestartIntvl",3,num_args);

  LuaCheckNilValue("chiLBSGroupsetSetGMRESRestartIntvl",L,1);
  LuaCheckNilValue("chiLBSGroupsetSetGMRESRestartIntvl",L,2);
  LuaCheckNilValue("chiLBSGroupsetSetGMRESRestartIntvl",L,3);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int restart_intvl = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiLBSGroupsetSetGMRESRestartIntvl";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiLBSGroupsetSetGMRESRestartIntvl";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  LBS_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiLBSGroupsetSetGMRESRestartIntvl";
    exit(EXIT_FAILURE);
  }

  //============================================= Bounds checking
  if (restart_intvl < 2)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid GMRES restart interval specified "
      << "in call to chiLBSGroupsetSetGMRESRestartIntvl. Must be >= 3.";
  }

  groupset->gmres_restart_intvl = restart_intvl;

  chi_log.Log(LOG_0)
    << "Groupset " << grpset_index << " GMRES restart interval set to "
    << "set to " << restart_intvl;

  return 0;
}

//###################################################################
/**Sets the Within-Group Diffusion Synthetic Acceleration parameters
 * for this groupset. If this call is being made then it is assumed
 * WGDSA is being applied.
 *
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param MaxIters int Maximum amount of iterations to use for WGDSA solvers.
                    Default 30.
\param ResTol float Residual tolerance to use for the WGDSA solve.

\param Verbose bool Optional flag indicating verbose output of WGDSA.
                    Default false.
\param PETSCString char Optional. Options string to be inserted
                        during initialization.



##_

Example:
\code
petsc_options =                  " -pc_hypre_boomeramg_strong_threshold 0.8"
petsc_options = petsc_options .. " -pc_hypre_boomeramg_max_levels 25"
chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false,petsc_options)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiLBSGroupsetSetWGDSA(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args < 4)
    LuaPostArgAmountError("chiLBSGroupsetSetWGDSA",4,num_args);

  LuaCheckNilValue("chiLBSGroupsetSetWGDSA",L,1);
  LuaCheckNilValue("chiLBSGroupsetSetWGDSA",L,2);
  LuaCheckNilValue("chiLBSGroupsetSetWGDSA",L,3);
  LuaCheckNilValue("chiLBSGroupsetSetWGDSA",L,4);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int max_iters = lua_tonumber(L,3);
  double resid_tol = lua_tonumber(L,4);
  bool verbose = false;
  const char* petsc_string = "";

  if (num_args >= 5)
    verbose = lua_toboolean(L,5);

  if (num_args == 6)
    petsc_string = lua_tostring(L,6);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiLBSGroupsetSetWGDSA";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiLBSGroupsetSetWGDSA";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  LBS_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiLBSGroupsetSetWGDSA";
    exit(EXIT_FAILURE);
  }

  groupset->apply_wgdsa     = true;
  groupset->wgdsa_max_iters = max_iters;
  groupset->wgdsa_tol       = resid_tol;
  groupset->wgdsa_verbose   = verbose;
  groupset->wgdsa_string    = std::string(petsc_string);

  chi_log.Log(LOG_0)
    << "Groupset " << grpset_index << " set to apply WGDSA with "
    << max_iters << " maximum iterations and a tolerance of "
    << resid_tol
    << ". PETSc-string: " << std::string(petsc_string);

  return 0;
}

//###################################################################
/**Sets the Two-Grid Diffusion Synthetic Acceleration parameters
 * for this groupset. If this call is being made then it is assumed
 * TGDSA is being applied.
 *
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param MaxIters int Maximum amount of iterations to use for TGDSA solvers.
                    Default 30.
\param ResTol float Residual tolerance to use for the TGDSA solve.

\param Verbose bool Optional flag indicating verbose output of TGDSA.
                    Default false.
\param PETSCString char Optional. Options string to be inserted
                        during initialization.



##_

Example:
\code
petsc_options =                  " -pc_hypre_boomeramg_strong_threshold 0.8"
petsc_options = petsc_options .. " -pc_hypre_boomeramg_max_levels 25"
chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false,petsc_options)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiLBSGroupsetSetTGDSA(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args < 4)
    LuaPostArgAmountError("TGDSA",4,num_args);

  LuaCheckNilValue("chiLBSGroupsetSetTGDSA",L,1);
  LuaCheckNilValue("chiLBSGroupsetSetTGDSA",L,2);
  LuaCheckNilValue("chiLBSGroupsetSetTGDSA",L,3);
  LuaCheckNilValue("chiLBSGroupsetSetTGDSA",L,4);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int max_iters = lua_tonumber(L,3);
  double resid_tol = lua_tonumber(L,4);
  bool verbose = false;
  const char* petsc_string = "";

  if (num_args >= 5)
    verbose = lua_toboolean(L,5);

  if (num_args == 6)
    petsc_string = lua_tostring(L,6);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiLBSGroupsetSetTGDSA";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiLBSGroupsetSetTGDSA";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  LBS_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiLBSGroupsetSetTGDSA";
    exit(EXIT_FAILURE);
  }

  groupset->apply_tgdsa     = true;
  groupset->tgdsa_max_iters = max_iters;
  groupset->tgdsa_tol       = resid_tol;
  groupset->tgdsa_verbose   = verbose;
  groupset->tgdsa_string    = std::string(petsc_string);

  chi_log.Log(LOG_0)
    << "Groupset " << grpset_index << " set to apply TGDSA with "
    << max_iters << " maximum iterations and a tolerance of "
    << resid_tol
    << ". PETSc-string: " << std::string(petsc_string);

  return 0;
}