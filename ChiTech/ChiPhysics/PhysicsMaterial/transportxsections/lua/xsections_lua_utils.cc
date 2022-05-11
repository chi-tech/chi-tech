#include "ChiLua/chi_lua.h"
#include<iostream>

#include "chi_runtime.h"

#include "ChiPhysics/chi_physics_namespace.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Creates a stand-alone transport cross-section.
 *
 *


\code
xs_grph_clean = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSet(xs_grph_clean,PDT_XSFILE,"xs_graphite_pure_116.data")

chiPhysicsMaterialSetProperty(materials[2],
                                    TRANSPORT_XSECTIONS,
                                    EXISTING,
                                    xs_grph_clean)
\endcode
 * \return Returns a handle to the cross-section.
 *
\ingroup LuaTransportXSs
 */
int chiPhysicsTransportXSCreate(lua_State* L)
{
  auto xs = std::make_shared<chi_physics::TransportCrossSections>();

  chi::trnsprt_xs_stack.push_back(xs);

  const size_t index = chi::trnsprt_xs_stack.size()-1;

  lua_pushinteger(L,static_cast<lua_Integer>(index));
  return 1;
}

//###################################################################
/**Sets the properties of a transport cross-section.

 \param XS_handle int Handle to the cross-section to be modified.
 \param OperationIndex int Method used for setting the xs property.
 \param Information varying Varying information depending on the operation.

 ##_

###OperationIndex
SINGLE_VALUE\n
Sets the property based on a single value. Requires a single value as additional
information. As a simple example consider the case where the user would like
to set a single constant thermal conductivity. This can be achieved with \n
FROM_ARRAY\n
Sets a property based on a Lua array indexed from 1 to N. Internally
will be converted to 0 to N-1. This method can be used to set mutligroup
cross-sections or sources.
\n
SIMPLEXS0\n
Makes a simple material with no transfer matrix just \f$\sigma_t \f$. Expects two
values: \n
 - int number of groups \f$G \f$,
 - float \f$\sigma_t \f$.

####_

SIMPLEXS1\n
Makes a simple material with isotropic transfer matrix (L=0)
and mostly down scattering but with a few of the last groups
subject to up-scattering. Expects three values
values: \n
 - int number of groups (\f$G \f$),
 - float \f$\sigma_t \f$,
 - float scattering to total ratio (\f$c \f$)

####_

PDT_XSFILE\n
Loads transport cross-sections from PDT type cross-section files. Expects
to be followed by a filepath specifying the xs-file. By default this routine
will attempt to build a transfer matrix from reaction type MT 2501, however,
an additional text field can be supplied specifying the transfer matrix to
 use.

####_

CHI_XSFILE\n
Loads transport cross-sections from CHI type cross-section files. Expects
to be followed by a filepath specifying the xs-file.


##_
### Example
Example lua code:
\code
graphite = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSet(graphite,"xs_3_170.data","2518")
\endcode
 *
 *
\ingroup LuaTransportXSs
 * \return */
int chiPhysicsTransportXSSet(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args < 3)
  {
    LuaPostArgAmountError("chiPhysicsTransportXSSet",3,num_args);
    exit(EXIT_FAILURE);
  }

  LuaCheckNilValue("chiPhysicsTransportXSSet",L,1);
  LuaCheckNilValue("chiPhysicsTransportXSSet",L,2);

  int handle = lua_tonumber(L,1);
  int operation_index = lua_tonumber(L,2);

  std::shared_ptr<chi_physics::TransportCrossSections> xs;
  try {
    xs = chi::GetStackItemPtr(chi::trnsprt_xs_stack, handle);
  }
  catch(const std::out_of_range& o){
    chi_log.Log(LOG_ALLERROR)
      << "ERROR: Invalid cross-section handle"
      << " in call to chiPhysicsTransportXSSet."
      << std::endl;
    exit(EXIT_FAILURE);
  }

  //========================== Process operation
  using OpType = chi_physics::OperationType;
  if (operation_index == static_cast<int>(OpType::SIMPLEXS0))
  {
    if (num_args!=4)
      LuaPostArgAmountError("chiPhysicsTransportXSSet",4,num_args);

    int    G       = lua_tonumber(L,3);
    double sigma_t = lua_tonumber(L,4);

    xs->MakeSimple0(G,sigma_t);
  }
  else if (operation_index == static_cast<int>(OpType::SIMPLEXS1))
  {
    if (num_args!=5)
      LuaPostArgAmountError("chiPhysicsTransportXSSet",5,num_args);

    int    G       = lua_tonumber(L,3);
    double sigma_t = lua_tonumber(L,4);
    double c       = lua_tonumber(L,5);

    xs->MakeSimple1(G,sigma_t,c);
  }
  else if (operation_index == static_cast<int>(OpType::PDT_XSFILE))
  {
    if (!((num_args>=3) && (num_args<=4)))
      LuaPostArgAmountError("chiPhysicsTransportXSSet",3,num_args);

    const char* file_name_c = lua_tostring(L,3);
    std::string MT_TRANSFER("2501");

    if (num_args == 4)
      MT_TRANSFER = std::string(lua_tostring(L,4));

    xs->MakeFromPDTxsFile(std::string(file_name_c),MT_TRANSFER);
  }
  else if (operation_index == static_cast<int>(OpType::CHI_XSFILE))
  {
    if (num_args != 3)
      LuaPostArgAmountError("chiPhysicsTransportXSSet",3,num_args);

    const char* file_name_c = lua_tostring(L,3);

    xs->MakeFromCHIxsFile(std::string(file_name_c));
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "Unsupported operation in "
      << "chiPhysicsTransportXSSet. " << operation_index
      << std::endl;
    exit(EXIT_FAILURE);
  }
  return 0;
}

//###################################################################
/**Obtains a lua table of all the cross-section values.

 \param XS_handle int Handle to the cross-section to be modified.

 ## _

 To print the contents of the table, execute the following:
\code
xs = chiPhysicsTransportXSGet(xs_handle)
for i,v in pairs(xs) do
    print(i,v)
end
\endcode
\ingroup LuaTransportXSs*/
int chiPhysicsTransportXSGet(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args < 1)
  {
    LuaPostArgAmountError(__FUNCTION__,1,num_args);
    exit(EXIT_FAILURE);
  }

  LuaCheckNilValue(__FUNCTION__,L,1);

  int handle = lua_tonumber(L,1);

  std::shared_ptr<chi_physics::TransportCrossSections> xs;
  try {
    xs = chi::GetStackItemPtr(chi::trnsprt_xs_stack, handle);
  }
  catch(const std::out_of_range& o){
    chi_log.Log(LOG_ALLERROR)
      << "ERROR: Invalid cross-section handle"
      << " in call to " << __FUNCTION__ << "."
      << std::endl;
    exit(EXIT_FAILURE);
  }

  xs->PushLuaTable(L);

  return 1;
}

//###################################################################
/**Makes a combined cross-section from multiple other cross-sections.

 \param Combinations table A lua-table with each element another table
                           containing a handle to an existing xs and a
                           scalar multiplier.

 ## _

###Example:
Example lua code:
\code
xs_1 = chiPhysicsTransportXSCreate()
xs_2 = chiPhysicsTransportXSCreate()
xs_3 = chiPhysicsTransportXSCreate()

chiPhysicsTransportXSSet(xs_1,PDT_XSFILE,"CHI_TEST/xs_graphite_pure.data")
chiPhysicsTransportXSSet(xs_2,PDT_XSFILE,"CHI_TEST/xs_3_170.data")
chiPhysicsTransportXSSet(xs_3,PDT_XSFILE,"CHI_TEST/xs_air50RH.data")

combo ={{xs_1, 0.5e5},
        {xs_2, 0.4e3},
        {xs_3, 0.3e2}}
aerated_graphite = chiPhysicsTransportXSMakeCombined(combo)


chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              EXISTING,aerated_graphite)
\endcode

 \return Returns a handle to another cross-section object that contains the
         desired combination.

\ingroup LuaTransportXSs
 */
int chiPhysicsTransportXSMakeCombined(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiPhysicsMakeCombinedTransportXS",1,num_args);

  if (!lua_istable(L,1))
  {
    chi_log.Log(LOG_ALLERROR)
      << "In call to chiPhysicsMakeCombinedTransportXS: "
      << "Argument must be a lua table.";
    exit(EXIT_FAILURE);
  }

  size_t table_len = lua_rawlen(L,1);

  std::vector<std::pair<int,double>> combinations;
  combinations.reserve(table_len);

  //======================================== Process table
  for (int v=0; v<table_len; ++v)
  {
    lua_pushnumber(L,v+1);
    lua_gettable(L,1);

    if (!lua_istable(L,-1))
    {
      chi_log.Log(LOG_ALLERROR)
        << "In call to chiPhysicsMakeCombinedTransportXS: "
        << "The elements of the supplied table must themselves also"
           "be lua tables of the xs handle and its scalar multiplier.";
      exit(EXIT_FAILURE);
    }

    lua_pushinteger(L,1);
    lua_gettable(L,-2);
    LuaCheckNilValue("chiPhysicsMakeCombinedTransportXS:A1:E1",L,-1);

    int handle = lua_tonumber(L,-1); lua_pop(L,1);

    lua_pushinteger(L,2);
    lua_gettable(L,-2);
    LuaCheckNilValue("chiPhysicsMakeCombinedTransportXS:A1:E2",L,-1);

    double scalar = lua_tonumber(L,-1); lua_pop(L,1);

    combinations.emplace_back(handle,scalar);
    lua_pop(L,1); //pop off table
  }

  //======================================== Print out table
  chi_log.Log(LOG_0) << "Generating XS with following combination:";
  for (auto& elem : combinations)
    chi_log.Log(LOG_0) << "  Element handle: " << elem.first
                       << " scalar value: " << elem.second;

  //======================================== Make the new cross-section
  auto new_xs = std::make_shared<chi_physics::TransportCrossSections>();

  new_xs->MakeCombined(combinations);

  chi::trnsprt_xs_stack.push_back(new_xs);
  lua_pushinteger(L,
      static_cast<lua_Integer>(chi::trnsprt_xs_stack.size())-1);

  return 1;
}

//###################################################################
/**Sets a combined cross-section from multiple other cross-sections. This
 function can be called multiple times on the same cross-section handle.

 \param XS_handle int Handle to the cross-section to be modified.
 \param Combinations table A lua-table with each element another table
                           containing a handle to an existing xs and a
                           scalar multiplier.

 ## _

###Example:
Example lua code:
\code
xs_1 = chiPhysicsTransportXSCreate()
xs_2 = chiPhysicsTransportXSCreate()
xs_3 = chiPhysicsTransportXSCreate()

chiPhysicsTransportXSSet(xs_1,PDT_XSFILE,"CHI_TEST/xs_graphite_pure.data")
chiPhysicsTransportXSSet(xs_2,PDT_XSFILE,"CHI_TEST/xs_3_170.data")
chiPhysicsTransportXSSet(xs_3,PDT_XSFILE,"CHI_TEST/xs_air50RH.data")

combo ={{xs_1, 0.5e5},
        {xs_2, 0.4e3},
        {xs_3, 0.3e2}}
xs_4 = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSetCombined(xs_4,combo)

\endcode
 *
\ingroup LuaTransportXSs
 * */
int chiPhysicsTransportXSSetCombined(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args < 2)
  {
    LuaPostArgAmountError(__FUNCTION__,2,num_args);
    exit(EXIT_FAILURE);
  }

  LuaCheckNilValue(__FUNCTION__,L,1);
  LuaCheckNilValue(__FUNCTION__,L,2);
  LuaCheckTableValue(__FUNCTION__ ,L,2);

  //======================================== Process handle
  int xs_handle = lua_tonumber(L,1);

  std::shared_ptr<chi_physics::TransportCrossSections> xs;
  try {
    xs = chi::GetStackItemPtr(chi::trnsprt_xs_stack, xs_handle);
  }
  catch(const std::out_of_range& o){
    chi_log.Log(LOG_ALLERROR)
      << "ERROR: Invalid cross-section handle"
      << " in call to " << __FUNCTION__ << "."
      << std::endl;
    exit(EXIT_FAILURE);
  }

  //======================================== Process table
  size_t table_len = lua_rawlen(L,2);

  std::vector<std::pair<int,double>> combinations;
  combinations.reserve(table_len);

  for (int v=0; v<table_len; ++v)
  {
    lua_pushnumber(L,v+1);
    lua_gettable(L,1);

    if (!lua_istable(L,-1))
    {
      chi_log.Log(LOG_ALLERROR)
        << "In call to " << __FUNCTION__ << ": "
        << "The elements of the supplied table must themselves also"
           "be lua tables of the xs handle and its scalar multiplier.";
      exit(EXIT_FAILURE);
    }

    lua_pushinteger(L,1);
    lua_gettable(L,-2);
    LuaCheckNilValue((std::string(__FUNCTION__) + ":A1:E1").c_str(),L,-1);

    int handle = lua_tonumber(L,-1); lua_pop(L,1);

    lua_pushinteger(L,2);
    lua_gettable(L,-2);
    LuaCheckNilValue((std::string(__FUNCTION__) + ":A1:E2").c_str(),L,-1);

    double scalar = lua_tonumber(L,-1); lua_pop(L,1);

    combinations.emplace_back(handle,scalar);
    lua_pop(L,1); //pop off table
  }

  //======================================== Print out table
  chi_log.Log(LOG_0) << "Setting XS with following combination:";
  for (auto& elem : combinations)
    chi_log.Log(LOG_0) << "  Element handle: " << elem.first
                       << " scalar value: " << elem.second;

  xs->MakeCombined(combinations);

  return 0;
}

//###################################################################
/** Exports a cross section to ChiTech format.
 *
\param XS_handle int Handle to the cross-section to be exported.
\param file_name string The name of the file to which the XS is to be exported.

\ingroup LuaTransportXSs
 */
int chiPhysicsTransportXSExportToChiTechFormat(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args != 2)
    LuaPostArgAmountError(__FUNCTION__,2,num_args);

  LuaCheckNilValue(__FUNCTION__,L,1);
  LuaCheckNilValue(__FUNCTION__,L,2);

  //======================================== Process handle
  int handle = lua_tonumber(L,1);

  std::shared_ptr<chi_physics::TransportCrossSections> xs;
  try {
    xs = chi::GetStackItemPtr(chi::trnsprt_xs_stack, handle);
  }
  catch(const std::out_of_range& o){
    chi_log.Log(LOG_ALLERROR)
      << "ERROR: Invalid cross-section handle"
      << " in call to " << __FUNCTION__ << "."
      << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string file_name = lua_tostring(L,2);

  xs->ExportToChiFormat(file_name);

  return 0;
}