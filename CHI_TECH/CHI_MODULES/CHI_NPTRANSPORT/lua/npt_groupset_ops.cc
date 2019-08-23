#include "../../../CHI_LUA/chi_lua.h"

#include "../chi_nptransport.h"
#include "../../../CHI_PHYSICS/chi_physics.h"
#include "../../../CHI_MATH/chi_math.h"
#include <chi_log.h>

extern CHI_PHYSICS chi_physics_handler;
extern CHI_MATH    chi_math_handler;
extern CHI_LOG     chi_log;

/** \defgroup LuaLBSGroupsets LBS Groupsets

The code below is an example of a complete specification of a groupset.

\code
--===================================== Setup physics
phys1 = chiNPTransportCreateSolver()
chiSolverAddRegion(phys1,region1)

chiNPTSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
chiNPTSetProperty(phys1,SCATTERING_ORDER,1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiNPTCreateGroup(phys1)
end

--========== ProdQuad
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)
pquad1 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,8, 8)

--========== Groupset def
gs0 = chiNPTCreateGroupset(phys1)

cur_gs = gs0
chiNPTGroupsetAddGroups(phys1,cur_gs,0,15)
chiNPTGroupsetSetQuadrature(phys1,cur_gs,pquad0)
chiNPTGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiNPTGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiNPTGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
chiNPTGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-4)
chiNPTGroupsetSetMaxIterations(phys1,cur_gs,300)
chiNPTGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)
chiNPTGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
chiNPTGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")
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
gs0 = chiNPTCreateGroupset(phys1)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiNPTCreateGroupset(lua_State *L)
{
  int solver_index = lua_tonumber(L,1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  CHI_NPTRANSPORT* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
    {
      solver = (CHI_NPTRANSPORT*)(psolver);
    }
    else
    {
      fprintf(stderr,"ERROR: Incorrect solver-type"
                     "in chiNPTCreateGroupset\n");
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiNPTCreateGroupset\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Create groupset
  NPT_GROUPSET* newgs = new NPT_GROUPSET;
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
grp[g] = chiNPTCreateGroup(phys1)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiNPTCreateGroup(lua_State *L)
{
  int solver_index = lua_tonumber(L,1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  CHI_NPTRANSPORT* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
    {
      solver = (CHI_NPTRANSPORT*)(psolver);
    }
    else
    {
      fprintf(stderr,"ERROR: Incorrect solver-type"
                     "in chiNPTCreateGroup\n");
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiNPTCreateGroup\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Create groupset
  NPT_GROUP* newgs = new NPT_GROUP;
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
    grp[g] = chiNPTCreateGroup(phys1)
end

chiNPTGroupsetAddGroups(phys1,cur_gs,0,15)
\endcode


\ingroup LuaLBSGroupsets
*/
int chiNPTGroupsetAddGroups(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 4)
  {
    fprintf(stderr,"ERROR: Invalid amount of arguments, %d,"
                   "in chiNPTGroupsetAddGroups"
                   " (4 required)\n",num_args);
    exit(EXIT_FAILURE);
  }
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int from = lua_tonumber(L,3);
  int to   = lua_tonumber(L,4);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  CHI_NPTRANSPORT* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
    {
      solver = (CHI_NPTRANSPORT*)(psolver);
    }
    else
    {
      fprintf(stderr,"ERROR: Incorrect solver-type in "
                     "chiNPTGroupsetAddGroups\n");
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiNPTGroupsetAddGroups\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  NPT_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to groupset"
                   "in chiNPTGroupsetAddGroups\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Add the groups
  if (to<from)
  {
    chi_log.Log(LOG_0ERROR)
    << "No groups added to groupset in chiNPTGroupsetAddGroups. "
       "This is triggered when groups are added with the \"to\" "
       "field being less than the \"from\" field.";
    exit(EXIT_FAILURE);
  }


  for (unsigned k=from; k<=to; k++)
  {
    NPT_GROUP* group;
    //================================= Check valid group
    try {
      group = solver->groups.at(k);
    }
    catch (std::out_of_range o)
    {
      fprintf(stderr,"ERROR: Invalid group added to groupset"
                     "in chiNPTGroupsetAddGroups\n");
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

chiNPTGroupsetSetQuadrature(phys1,cur_gs,pquad0)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiNPTGroupsetSetQuadrature(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiNPTGroupsetSetQuadrature",3,num_args);

  LuaCheckNilValue("chiNPTGroupsetSetQuadrature",L,1);
  LuaCheckNilValue("chiNPTGroupsetSetQuadrature",L,2);
  LuaCheckNilValue("chiNPTGroupsetSetQuadrature",L,3);


  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int prquad_index = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  CHI_NPTRANSPORT* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
    {
      solver = (CHI_NPTRANSPORT*)(psolver);
    }
    else
    {
      fprintf(stderr,"ERROR: Incorrect solver-type"
                     "in chiNPTGroupsetSetQuadrature\n");
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiNPTGroupsetSetQuadrature\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  NPT_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to groupset"
                   "in chiNPTGroupsetSetQuadrature\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to quadrature
  CHI_PRODUCT_QUADRATURE* prodquad;
  try{
    prodquad = chi_math_handler.product_quadratures.at(prquad_index);
  }
  catch (std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to Product Quadrature"
                   "in chiNPTGroupsetSetQuadrature\n");
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
chiNPTGroupsetSetAngleAggDiv(phys1,cur_gs,1)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiNPTGroupsetSetAngleAggDiv(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiNPTGroupsetSetAngleAggDiv",3,num_args);

  LuaCheckNilValue("chiNPTGroupsetSetAngleAggDiv",L,1);
  LuaCheckNilValue("chiNPTGroupsetSetAngleAggDiv",L,2);
  LuaCheckNilValue("chiNPTGroupsetSetAngleAggDiv",L,3);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int num_div = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  CHI_NPTRANSPORT* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
    {
      solver = (CHI_NPTRANSPORT*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiNPTGroupsetSetAngleAggDiv";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiNPTGroupsetSetAngleAggDiv";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  NPT_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiNPTGroupsetSetAngleAggDiv";
    exit(EXIT_FAILURE);
  }

  //============================================= Bounds checking
  if (num_div <= 0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid number of divisions "
      << "in call to chiNPTGroupsetSetAngleAggDiv. Must be >= 1.";
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
chiNPTGroupsetSetGroupSubsets(phys1,cur_gs,1)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiNPTGroupsetSetGroupSubsets(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiNPTGroupsetSetGroupSubsets",3,num_args);

  LuaCheckNilValue("chiNPTGroupsetSetGroupSubsets",L,1);
  LuaCheckNilValue("chiNPTGroupsetSetGroupSubsets",L,2);
  LuaCheckNilValue("chiNPTGroupsetSetGroupSubsets",L,3);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int num_div = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  CHI_NPTRANSPORT* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
    {
      solver = (CHI_NPTRANSPORT*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiNPTGroupsetSetGroupSubsets";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiNPTGroupsetSetGroupSubsets";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  NPT_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiNPTGroupsetSetGroupSubsets";
    exit(EXIT_FAILURE);
  }

  //============================================= Bounds checking
  if (num_div <= 0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid number of subsets "
      << "in call to chiNPTGroupsetSetGroupSubsets. Must be >= 1.";
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
chiNPTGroupsetSetIterativeMethod(phys1,cur_gs,NPT_CLASSICRICHARDSON)
chiNPTGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiNPTGroupsetSetIterativeMethod(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiNPTGroupsetSetIterativeMethod",3,num_args);

  LuaCheckNilValue("chiNPTGroupsetSetIterativeMethod",L,1);
  LuaCheckNilValue("chiNPTGroupsetSetIterativeMethod",L,2);
  LuaCheckNilValue("chiNPTGroupsetSetIterativeMethod",L,3);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int iter_method  = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  CHI_NPTRANSPORT* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
    {
      solver = (CHI_NPTRANSPORT*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiNPTGroupsetSetGroupSubsets";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiNPTGroupsetSetGroupSubsets";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  NPT_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiNPTGroupsetSetGroupSubsets";
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
      << "chiNPTGroupsetSetIterativeMethod.";
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
chiNPTGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-4)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiNPTGroupsetSetResidualTolerance(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiNPTGroupsetSetResidualTolerance",3,num_args);

  LuaCheckNilValue("chiNPTGroupsetSetResidualTolerance",L,1);
  LuaCheckNilValue("chiNPTGroupsetSetResidualTolerance",L,2);
  LuaCheckNilValue("chiNPTGroupsetSetResidualTolerance",L,3);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  double resid_tol = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  CHI_NPTRANSPORT* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
    {
      solver = (CHI_NPTRANSPORT*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiNPTGroupsetSetGroupSubsets";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiNPTGroupsetSetGroupSubsets";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  NPT_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiNPTGroupsetSetGroupSubsets";
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
chiNPTGroupsetSetMaxIterations(phys1,cur_gs,200)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiNPTGroupsetSetMaxIterations(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiNPTGroupsetSetMaxIterations",3,num_args);

  LuaCheckNilValue("chiNPTGroupsetSetMaxIterations",L,1);
  LuaCheckNilValue("chiNPTGroupsetSetMaxIterations",L,2);
  LuaCheckNilValue("chiNPTGroupsetSetMaxIterations",L,3);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int num_iter = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  CHI_NPTRANSPORT* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
    {
      solver = (CHI_NPTRANSPORT*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiNPTGroupsetSetMaxIterations";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiNPTGroupsetSetMaxIterations";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  NPT_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiNPTGroupsetSetMaxIterations";
    exit(EXIT_FAILURE);
  }

  //============================================= Bounds checking
  if (num_iter < 0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid number of iterations "
      << "in call to chiNPTGroupsetSetMaxIterations. Must be >= 0.";
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
chiNPTGroupsetSetGMRESRestartIntvl(phys1,cur_gs,15)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiNPTGroupsetSetGMRESRestartIntvl(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiNPTGroupsetSetGMRESRestartIntvl",3,num_args);

  LuaCheckNilValue("chiNPTGroupsetSetGMRESRestartIntvl",L,1);
  LuaCheckNilValue("chiNPTGroupsetSetGMRESRestartIntvl",L,2);
  LuaCheckNilValue("chiNPTGroupsetSetGMRESRestartIntvl",L,3);
  int solver_index = lua_tonumber(L,1);
  int grpset_index = lua_tonumber(L,2);
  int restart_intvl = lua_tonumber(L,3);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  CHI_NPTRANSPORT* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
    {
      solver = (CHI_NPTRANSPORT*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiNPTGroupsetSetGMRESRestartIntvl";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiNPTGroupsetSetGMRESRestartIntvl";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  NPT_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiNPTGroupsetSetGMRESRestartIntvl";
    exit(EXIT_FAILURE);
  }

  //============================================= Bounds checking
  if (restart_intvl < 2)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid GMRES restart interval specified "
      << "in call to chiNPTGroupsetSetGMRESRestartIntvl. Must be >= 3.";
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
chiNPTGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false,petsc_options)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiNPTGroupsetSetWGDSA(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args < 4)
    LuaPostArgAmountError("chiNPTGroupsetSetWGDSA",4,num_args);

  LuaCheckNilValue("chiNPTGroupsetSetWGDSA",L,1);
  LuaCheckNilValue("chiNPTGroupsetSetWGDSA",L,2);
  LuaCheckNilValue("chiNPTGroupsetSetWGDSA",L,3);
  LuaCheckNilValue("chiNPTGroupsetSetWGDSA",L,4);
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
  CHI_NPTRANSPORT* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
    {
      solver = (CHI_NPTRANSPORT*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiNPTGroupsetSetWGDSA";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiNPTGroupsetSetWGDSA";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  NPT_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiNPTGroupsetSetWGDSA";
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
chiNPTGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false,petsc_options)
\endcode

\ingroup LuaLBSGroupsets
*/
int chiNPTGroupsetSetTGDSA(lua_State *L)
{
  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args < 4)
    LuaPostArgAmountError("TGDSA",4,num_args);

  LuaCheckNilValue("chiNPTGroupsetSetTGDSA",L,1);
  LuaCheckNilValue("chiNPTGroupsetSetTGDSA",L,2);
  LuaCheckNilValue("chiNPTGroupsetSetTGDSA",L,3);
  LuaCheckNilValue("chiNPTGroupsetSetTGDSA",L,4);
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
  CHI_NPTRANSPORT* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
    {
      solver = (CHI_NPTRANSPORT*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Incorrect solver-type "
        << "in call to chiNPTGroupsetSetTGDSA";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to solver "
      << "in call to chiNPTGroupsetSetTGDSA";
    exit(EXIT_FAILURE);
  }

  //============================================= Obtain pointer to groupset
  NPT_GROUPSET* groupset;
  try{
    groupset = solver->group_sets.at(grpset_index);
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to groupset "
      << "in call to chiNPTGroupsetSetTGDSA";
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