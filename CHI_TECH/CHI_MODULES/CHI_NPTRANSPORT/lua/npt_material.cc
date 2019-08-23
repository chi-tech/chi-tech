#include "../../../CHI_LUA/chi_lua.h"

#include "../chi_nptransport.h"
#include "../../../CHI_PHYSICS/chi_physics.h"

extern CHI_PHYSICS chi_physics_handler;


////###################################################################
///**Create an empty material.
//\param SolverIndex int Handle to the solver for which the material
// is to be created.
//\ingroup LuaNPT
//*/
//int chiNPTCreateMaterial(lua_State *L)
//{
//  int solver_index = lua_tonumber(L,1);
//  char* matnam = (char*)lua_tostring(L,2);
//
//  //============================================= Get pointer to solver
//  chi_physics::Solver* psolver;
//  CHI_NPTRANSPORT* solver;
//  try{
//    psolver = chi_physics_handler.solver_stack.at(solver_index);
//
//    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
//    {
//      solver = (CHI_NPTRANSPORT*)(psolver);
//    }
//    else
//    {
//      fprintf(stderr,"ERROR: Incorrect solver-type"
//                     "in chiNPTCreateGroupset\n");
//      exit(EXIT_FAILURE);
//    }
//  }
//  catch(std::out_of_range o)
//  {
//    fprintf(stderr,"ERROR: Invalid handle to solver"
//                   "in chiNPTCreateGroupset\n");
//    exit(EXIT_FAILURE);
//  }
//
//  //============================================= Create material
//  NPT_MATERIAL* newmat = new NPT_MATERIAL;
//  newmat->name = std::string(matnam);
//  solver->materials.push_back(newmat);
//  newmat->id = solver->materials.size()-1;
//
//  lua_pushnumber(L,newmat->id);
//
//  return 1;
//}
//
//
//
////###################################################################
///**Set the material as a pure absorber.
//\param SolverIndex int Handle to the solver for which the material
// is to be created.
//\param MaterialIndex int Handle to the material in question.
//\param NumberOfGroups int Number of groups for the material.
//\param ScatOrder int Maximum scattering order of the material.
//\param Sigma_t=1.0 double (Optional) Total cross-section to be used.
//\ingroup LuaNPT
//*/
//int chiNPTMaterialSetPureAbsorber(lua_State *L)
//{
//  //============================================= Get arguments
//  int num_args = lua_gettop(L);
//  if (num_args<4)
//  {
//    fprintf(stderr,"ERROR: Invalid amount of arguments, %d,"
//                   "in chiNPTMaterialSetPureAbsorber"
//                   " (4 min required)\n",num_args);
//    exit(EXIT_FAILURE);
//  }
//  int solver_index = lua_tonumber(L,1);
//  int material_index = lua_tonumber(L,2);
//  int number_of_grps = lua_tonumber(L,3);
//  int scat_order     = lua_tonumber(L,4);
//  double sigma_t = 1.0;
//
//  if (num_args == 5)
//  {
//    sigma_t = lua_tonumber(L,5);
//  }
//
//  //============================================= Get pointer to solver
//  chi_physics::Solver* psolver;
//  CHI_NPTRANSPORT* solver;
//  try{
//    psolver = chi_physics_handler.solver_stack.at(solver_index);
//
//    if (typeid(*psolver) == typeid(CHI_NPTRANSPORT))
//    {
//      solver = (CHI_NPTRANSPORT*)(psolver);
//    }
//    else
//    {
//      fprintf(stderr,"ERROR: Incorrect solver-type"
//                     "in chiNPTMaterialSetPureAbsorber\n");
//      exit(EXIT_FAILURE);
//    }
//  }
//  catch(std::out_of_range o)
//  {
//    fprintf(stderr,"ERROR: Invalid handle to solver"
//                   "in chiNPTMaterialSetPureAbsorber\n");
//    exit(EXIT_FAILURE);
//  }
//
//  //============================================= Get pointer to material
//  NPT_MATERIAL* material;
//  try{
//    material = solver->materials.at(material_index);
//  }
//  catch (std::out_of_range o)
//  {
//    fprintf(stderr,"ERROR: Invalid handle to material"
//                   "in chiNPTMaterialSetPureAbsorber\n");
//    exit(EXIT_FAILURE);
//  }
//
//  //============================================= Set as pure absorber
//  material->SetAsPureAbsorber(number_of_grps,scat_order,sigma_t);
//
//  return 0;
//}