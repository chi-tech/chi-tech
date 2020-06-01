#include "../../ChiLua/chi_lua.h"
#include<iostream>
#include "../chi_physics.h"
#include "../PhysicsMaterial/chi_physicsmaterial.h"
#include "ChiPhysics/PhysicsMaterial/property01_scalarvalue.h"
#include "ChiPhysics/PhysicsMaterial/property10_transportxsections.h"
#include "ChiPhysics/PhysicsMaterial/property11_isotropic_mg_src.h"

extern ChiPhysics&  chi_physics_handler;

#include <chi_log.h>

extern ChiLog& chi_log;


//#############################################################################
/** Adds a material property to a material.
 *
\param MaterialHandle int Index to the reference material.
\param PropertyIndex int Property index.

##_

###PropertyIndex\n

SCALAR_VALUE\n
 Simple scalar value property.\n\n

TRANSPORT_XSECTIONS\n
 Multi-group transport cross-section supporting numerous features.\n\n

ISOTROPIC_MG_SOURCE\n
 Isotropic Multigroup Source.\n\n

### Developer Info
Checklist for adding a new material property:
 - Create your property class in its own header file. i.e.
   "ChiPhysics/PhysicsMaterial/property_xx_myproperty.h"
 - Add the property to the physics namespace
   ("ChiPhysics/chi_physics_namespace.h"). Make sure to derive from the base
   class.
 - Go define the integer to be associated with your new property in
   chi_physicsmaterial.h
 - Include the header file for your property in this file (i.e. at the top).
 - Add this property integer in the lua register (ChiLua/chi_lua_register.h).
   For testing you can just use the integer value but eventually you'd want
   to supply an easier way for users to enter it.
 - Add another else-if for your property. Just have a look at how the others
   were done, it should be intuitive enough.

##_

### Example\n
Example lua code:
\code
chiPhysicsMaterialAddProperty(materials[i],TRANSPORT_XSECTIONS)
\endcode

\ingroup LuaPhysicsMaterials
\author Jan*/
int chiPhysicsMaterialAddProperty(lua_State *L)
{
  int numArgs = lua_gettop(L);

  if (!((numArgs>=2) && (numArgs<=3)))
  {
    chi_log.Log(LOG_0ERROR) << "Incorrect amount of arguments "
                               "in chiPhysicsMaterialAddProperty";
    exit(EXIT_FAILURE);
  }

  int material_index = lua_tonumber(L,1);
  int property_index = lua_tonumber(L,2);

  const char* provided_name = "";
  if (numArgs == 3)
  {
    provided_name = lua_tostring(L,3);
  }

  //============================================= Get reference to material
  chi_physics::Material* cur_material;
  try {
    cur_material = chi_physics_handler.material_stack.at(material_index);
  }
  catch(const std::out_of_range& o){
    chi_log.Log(LOG_0ERROR) << "ERROR: Invalid material handle." << std::endl;
    exit(EXIT_FAILURE);
  }

  //============================================= Process property
  using MatProperty = chi_physics::PropertyType;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCALAR_VALUE
  if (property_index == static_cast<int>(MatProperty::SCALAR_VALUE))
  {
    //Duplicates are allowed

    auto prop = new chi_physics::ScalarValue;

    prop->property_name = std::string("Property ") +
                          std::to_string(cur_material->properties.size());

    if (numArgs == 3)
      prop->property_name = std::string(provided_name);

    cur_material->properties.push_back(prop);
    chi_log.Log(LOG_0VERBOSE_1) << "Scalar Value Property added to material"
                                 " at index " << material_index;
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSPORT_XSECTIONS
  else if (property_index == static_cast<int>(MatProperty::TRANSPORT_XSECTIONS))
  {
    //================================= Check for duplicate
    for (int p=0; p<cur_material->properties.size(); p++)
    {
      if (cur_material->properties[p]->Type() ==
            MatProperty::TRANSPORT_XSECTIONS)
      {
        chi_log.Log(LOG_0ERROR)    << "Material " << material_index << " \""
                                   << cur_material->name << "\""
                                   << " already has property "
                                      "TRANSPORT_XSECTIONS"
                                   << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    auto prop = new chi_physics::TransportCrossSections;

    prop->property_name = std::string("Property ") +
                          std::to_string(cur_material->properties.size());

    if (numArgs == 3)
      prop->property_name = std::string(provided_name);

    cur_material->properties.push_back(prop);
    chi_log.Log(LOG_0VERBOSE_1) << "Transport cross-sections added to material"
                                 " at index " << material_index;

    chi_physics_handler.trnsprt_xs_stack.push_back(prop);

    size_t index = chi_physics_handler.trnsprt_xs_stack.size()-1;

    lua_pushnumber(L,index);
    return 1;
  }
  else if (property_index == static_cast<int>(MatProperty::ISOTROPIC_MG_SOURCE))
  {
    //================================= Check for duplicate
    for (int p=0; p<cur_material->properties.size(); p++)
    {
      if (cur_material->properties[p]->Type() ==
            MatProperty::ISOTROPIC_MG_SOURCE)
      {
        chi_log.Log(LOG_0ERROR)    << "Material " << material_index << " \""
                                   << cur_material->name << "\""
                                   << " already has property "
                                      "ISOTROPIC_MG_SOURCE "
                                   << property_index
                                   << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    auto prop = new chi_physics::IsotropicMultiGrpSource;

    prop->property_name = std::string("Property ") +
                          std::to_string(cur_material->properties.size());

    if (numArgs == 3)
      prop->property_name = std::string(provided_name);

    cur_material->properties.push_back(prop);
    chi_log.Log(LOG_0VERBOSE_1) << "Isotropic Multigroup Source added to material"
                                 " at index " << material_index;
  }
  else
  {
    chi_log.Log(LOG_0ERROR)
      << "Unsupported property type in call to chiPhysicsMaterialAddProperty.";
    exit(EXIT_FAILURE);
  }


  return 0;
}


//#############################################################################
/** Sets a material property for a given material.
 *
\param MaterialHandle int Index to the reference material.
\param PropertyIndex int Property index. Or name of property.
\param OperationIndex int Method used for setting the material property.
\param Information varying Varying information depending on the operation.

##_

###PropertyIndex\n
SCALAR_VALUE         =  Basic Scalar value property.\n
TRANSPORT_XSECTIONS   =  Multi-group transport cross-section supporting numerous
                        features.\n
ISOTROPIC_MG_SOURCE = Isotropic Multigroup Source.\n

###OperationIndex\n
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
will attempt to build a transfer matrix from reaction type MT2501, however,
an additional text field can be supplied specifying the transfer matrix to
 use.

####_

CHI_XSFILE\n
Loads transport cross-sections from CHI type cross-section files. Expects
to be followed by a filepath specifying the xs-file. 

####_

EXISTING\n
Supply handle to an existing cross-section and simply swap them out.

\code
chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              PDT_XSFILE,
                              "xs_3_170.data",
                              "2518")
\endcode

##_

### Example 1
Simple temperature independent thermal conductivity:
\code
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
chiPhysicsMaterialAddProperty(materials[0],THERMAL_CONDUCTIVITY)
chiPhysicsMaterialSetProperty(materials[0],THERMAL_CONDUCTIVITY,SINGLE_VALUE,13.7)
\endcode

where the thermal conductivity has been set to 13.7.\n

### Example 2
Isotropic Multigroup source set from a lua table/array (12 groups):
\code
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)

num_groups = 12
src={}
for g=1,num_groups do
    src[g] = 0.0
end
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
\endcode

### Developer Info
Checklist for adding a new material property:
 - Make sure you followed the steps depicted in the developer info section for
   the ChiLua::chiPhysicsMaterialAddProperty function.
 - Now under the "If user supplied name then find property index"-section
   add the appropriate code for setting the property index.
 - Add an else-if block for your property similar to the others. It should be
   intuitive if you look at the others.
 - Remember to add the filtering section if you need to support multiple type
   properties.

\ingroup LuaPhysicsMaterials
\author Jan*/
int chiPhysicsMaterialSetProperty(lua_State *L)
{
  int numArgs = lua_gettop(L);

  if (numArgs<3)
  {
    chi_log.Log(LOG_0ERROR) << "Incorrect amount of arguments "
                               "in chiPhysicsMaterialSetProperty";
    exit(EXIT_FAILURE);
  }

  int material_index = lua_tonumber(L,1);
  int property_index = -1;
  std::string property_index_name;
  if (lua_isnumber(L,2))
    property_index = lua_tonumber(L,2);
  else
  {
    const char* temp_name = lua_tostring(L,2);
    property_index_name = std::string(temp_name);
  }

  int operation_index = lua_tonumber(L,3);

  //============================================= Get reference to material
  chi_physics::Material* cur_material;
  try {
    cur_material = chi_physics_handler.material_stack.at(material_index);
  }
  catch(const std::out_of_range& o){
    chi_log.Log(LOG_ALLERROR) << "ERROR: Invalid material handle." << std::endl;
    exit(EXIT_FAILURE);
  }

  //============================================= If user supplied name then
  //                                              find property index
  if (!lua_isnumber(L,2))
  {
    for (auto property : cur_material->properties)
      if (property->property_name == property_index_name)
        property_index = static_cast<int>(property->Type());
  }

  //============================================= Process property
  using MatProperty = chi_physics::PropertyType;
  using OpType = chi_physics::OperationType;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCALAR_VALUE
  if (property_index == static_cast<int>(MatProperty::SCALAR_VALUE))
  {
    int location_of_prop = -1;
    //================================= Check if the material has this property
    if (lua_isnumber(L,2))
    {
      for (int p=0; p<cur_material->properties.size(); p++)
        if (cur_material->properties[p]->Type() == MatProperty::SCALAR_VALUE)
          location_of_prop = p;
    }
    else
    {
      for (int p=0; p<cur_material->properties.size(); p++)
        if (cur_material->properties[p]->property_name == property_index_name)
          location_of_prop = p;
    }


    //================================= If the property is valid
    if (location_of_prop>=0)
    {
      auto prop =
        (chi_physics::ScalarValue*)cur_material->properties[location_of_prop];

      //========================== Process operation
      if (operation_index == static_cast<int>(OpType::SINGLE_VALUE))
      {
        double value = lua_tonumber(L,4);
        prop->value = value;
        chi_log.Log(LOG_0VERBOSE_1) << "Scalar value for material"
                                     " at index " << material_index
                                  << " set to " << value;
      }
      else
      {
        chi_log.Log(LOG_0ERROR) << "ERROR: Unsupported operation for "
                                   "SCALAR_VALUE." << std::endl;
        exit(EXIT_FAILURE);
      }

    }
    else
    {
      chi_log.Log(LOG_0ERROR) << "ERROR: Material has no property "
                                 "SCALAR_VALUE." << std::endl;
      exit(EXIT_FAILURE);
    }
  }//if scalar value
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSPORT_XSECTIONS
  else if (property_index == static_cast<int>(MatProperty::TRANSPORT_XSECTIONS))
  {
    int location_of_prop = -1;
    //================================= Check if the material has this property
    if (lua_isnumber(L,2))
    {
      for (int p=0; p<cur_material->properties.size(); p++)
      {
        if (cur_material->properties[p]->Type() ==
              MatProperty::TRANSPORT_XSECTIONS)
        {
          location_of_prop = p;
        }
      }
    }
    else
    {
      for (int p=0; p<cur_material->properties.size(); p++)
      {
        if (cur_material->properties[p]->property_name == property_index_name)
        {
          location_of_prop = p;
        }
      }
    }

    //================================= If the property is valid
    if (location_of_prop>=0)
    {
      auto prop = (chi_physics::TransportCrossSections*)
                  cur_material->properties[location_of_prop];

      //========================== Process operation
      if (operation_index == static_cast<int>(OpType::SIMPLEXS0))
      {
        if (numArgs!=5)
          LuaPostArgAmountError("chiPhysicsMaterialSetProperty",5,numArgs);

        int    G       = lua_tonumber(L,4);
        double sigma_t = lua_tonumber(L,5);

        prop->MakeSimple0(G,sigma_t);
      }
      else if (operation_index == static_cast<int>(OpType::SIMPLEXS1))
      {
        if (numArgs!=6)
          LuaPostArgAmountError("chiPhysicsMaterialSetProperty",6,numArgs);

        int    G       = lua_tonumber(L,4);
        double sigma_t = lua_tonumber(L,5);
        double c       = lua_tonumber(L,6);

        prop->MakeSimple1(G,sigma_t,c);
      }
      else if (operation_index == static_cast<int>(OpType::PDT_XSFILE))
      {
        if (!((numArgs>=4) && (numArgs<=5)))
          LuaPostArgAmountError("chiPhysicsMaterialSetProperty",4,numArgs);

        const char* file_name_c = lua_tostring(L,4);
        std::string MT_TRANSFER("2501");

        if (numArgs == 5)
          MT_TRANSFER = std::string(lua_tostring(L,5));

        prop->MakeFromPDTxsFile(std::string(file_name_c),MT_TRANSFER);
      }
      else if (operation_index == static_cast<int>(OpType::CHI_XSFILE))
      {
        if (numArgs != 4)
          LuaPostArgAmountError("chiPhysicsMaterialSetProperty",4,numArgs);

        const char* file_name_c = lua_tostring(L,4);

        prop->MakeFromCHIxsFile(std::string(file_name_c));
      }
      else if (operation_index == static_cast<int>(OpType::EXISTING))
      {
        if (numArgs != 4)
          LuaPostArgAmountError("chiPhysicsMaterialSetProperty",4,numArgs);

        int handle = lua_tonumber(L,4);

        chi_physics::TransportCrossSections* xs;
        try {
          xs = chi_physics_handler.trnsprt_xs_stack.at(handle);
        }
        catch(const std::out_of_range& o){
          chi_log.Log(LOG_ALLERROR)
            << "ERROR: Invalid cross-section handle"
            << " in call to chiPhysicsMaterialSetProperty."
            << std::endl;
          exit(EXIT_FAILURE);
        }
//        auto old_prop = prop;
        prop = xs;

        cur_material->properties[location_of_prop] = prop;

//        delete old_prop; //Still debating if this should be deleted
      }
      else
      {
        chi_log.Log(LOG_ALLERROR) << "Unsupported operation for "
                                   "TRANSPORT_XSECTIONS." << std::endl;
        exit(EXIT_FAILURE);
      }

    }
    else
    {
      chi_log.Log(LOG_ALLERROR) << "Material has no property "
                                 "TRANSPORT_XSECTIONS." << std::endl;
      exit(EXIT_FAILURE);
    }
  }//if thermal conductivity
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ISOTROPIC_MG_SOURCE
  else if (property_index == static_cast<int>(MatProperty::ISOTROPIC_MG_SOURCE))
  {
    int location_of_prop = -1;
    //================================= Check if the material has this property
    if (lua_isnumber(L,2))
    {
      for (int p=0; p<cur_material->properties.size(); p++)
      {
        if (cur_material->properties[p]->Type() ==
              MatProperty::ISOTROPIC_MG_SOURCE)
        {
          location_of_prop = p;
        }
      }
    }
    else
    {
      for (int p=0; p<cur_material->properties.size(); p++)
      {
        if (cur_material->properties[p]->property_name == property_index_name)
        {
          location_of_prop = p;
        }
      }
    }

    //================================= If the property is valid
    if (location_of_prop>=0)
    {
      auto prop =
        (chi_physics::IsotropicMultiGrpSource*)
          cur_material->properties[location_of_prop];


      if (operation_index == static_cast<int>(OpType::SINGLE_VALUE))
      {
        if (numArgs!=4)
          LuaPostArgAmountError("chiPhysicsMaterialSetProperty",4,numArgs);

        double value = lua_tonumber(L,4);

        prop->source_value_g.resize(1,value);
        chi_log.Log(LOG_0VERBOSE_1) << "Isotropic Multigroup Source value "
                                     "for material"
                                     " at index " << material_index
                                  << " set to " << value;
      }
      else if (operation_index == static_cast<int>(OpType::FROM_ARRAY))
      {
        if (numArgs!=4)
          LuaPostArgAmountError("chiPhysicsMaterialSetProperty",4,numArgs);

        if (!lua_istable(L,4))
        {
          chi_log.Log(LOG_ALLERROR)
            << "In call to chiPhysicsMaterialSetProperty: "
            << "Material \"" << cur_material->name
            << "\", when setting "
            << "ISOTROPIC_MG_SOURCE using operation FROM_ARRAY, the fourth "
               "argument was detected not to be a lua table.";
          exit(EXIT_FAILURE);
        }

        int table_len = lua_rawlen(L,4);

        std::vector<double> values(table_len,0.0);
        for (int g=0; g<table_len; g++)
        {
          lua_pushnumber(L,g+1);
          lua_gettable(L,4);
          values[g] = lua_tonumber(L,-1);
          lua_pop(L,1);
        }

        prop->source_value_g.resize(table_len,0.0);
        std::copy(values.begin(),values.end(),prop->source_value_g.begin());
        chi_log.Log(LOG_0VERBOSE_1) << "Isotropic Multigroup Source populated "
                                  << " with " << table_len << " values";
      }
      else
      {
        chi_log.Log(LOG_ALLERROR) << "Unsupported operation for "
                                     "ISOTROPIC_MG_SOURCE." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    else
    {
      chi_log.Log(LOG_ALLERROR) << "Material \"" << cur_material->name
                                << "\" has no property "
                                   "ISOTROPIC_MG_SOURCE." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    chi_log.Log(LOG_ALLERROR) << "Unsupported material property specified in "
                               "call to chiPhysicsMaterialSetProperty."
                               << property_index
                               << std::endl;
    exit(EXIT_FAILURE);
  }


  return 0;
}

//###################################################################
/** Returns a rich lua data-structure of the required property.
 *
\param MaterialHandle int Index to the reference material.
\param PropertyIndex int Property index. Or name of property.


\ingroup LuaPhysicsMaterials
\return Lua table of the desired property.

*/
int chiPhysicsMaterialGetProperty(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("chiPhysicsMaterialGetProperty",2,num_args);

  int material_index = lua_tonumber(L,1);
  int property_index = -1;
  std::string property_index_name;
  if (lua_isnumber(L,2))
    property_index = lua_tonumber(L,2);
  else
  {
    const char* temp_name = lua_tostring(L,2);
    property_index_name = std::string(temp_name);
  }

  //============================================= Get reference to material
  chi_physics::Material* cur_material;
  try {
    cur_material = chi_physics_handler.material_stack.at(material_index);
  }
  catch(const std::out_of_range& o){
    chi_log.Log(LOG_ALLERROR) << "ERROR: Invalid material handle." << std::endl;
    exit(EXIT_FAILURE);
  }

  //============================================= If user supplied name then
  //                                              find property index
  if (!lua_isnumber(L,2))
  {
    for (auto property : cur_material->properties)
      if (property->property_name == property_index_name)
        property_index = static_cast<int>(property->Type());
  }

  //============================================= Process property
  bool property_polulated = false;
  for (auto property : cur_material->properties)
  {
    if (static_cast<int>(property->Type()) == property_index)
    {
      property->PushLuaTable(L);
      property_polulated = true;
    }
  }


  if (not property_polulated)
  {
    chi_log.Log(LOG_ALLERROR) << "Invalid material property specified in "
                                 "call to chiPhysicsMaterialGetProperty."
                              << property_index
                              << std::endl;
    exit(EXIT_FAILURE);
  }

  return 1;
}

//###################################################################
/**Creates a stand-alone transport cross-section.
 *
 *
\ingroup LuaPhysicsMaterials

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
 * */
int chiPhysicsTransportXSCreate(lua_State* L)
{
  auto xs = new chi_physics::TransportCrossSections;

  chi_physics_handler.trnsprt_xs_stack.push_back(xs);

  size_t index = chi_physics_handler.trnsprt_xs_stack.size()-1;

  lua_pushnumber(L,index);
  return 1;
}

//###################################################################
/**Sets the properties of a transport cross-section.
 *
 * \param XS_handle int Handle to the cross-section to be modified.
 * \param OperationIndex int Method used for setting the xs property.
 * \param Information varying Varying information depending on the operation.
 *
 * ##_
 *
###OperationIndex\n
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

##_
### Example\n
Example lua code:
\code
graphite = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSet(graphite,"xs_3_170.data","2518")
\endcode
 *
 *
\ingroup LuaPhysicsMaterials
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

  chi_physics::TransportCrossSections* xs;
  try {
    xs = chi_physics_handler.trnsprt_xs_stack.at(handle);
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
/**Makes a combined cross-section from multiple other cross-sections.
 *
 * \param Combinations table A lua-table with each element another table
 *                           containing a handle to an existing xs and a
 *                           scalar multiplier.
 *
 * ## _
 *
###Example:\n
Example lua code:
\code
xs_1 = chiPhysicsTransportXSCreate()
xs_2 = chiPhysicsTransportXSCreate()
xs_3 = chiPhysicsTransportXSCreate()

chiPhysicsTransportXSSet(xs_1,PDT_XSFILE,"CHI_TEST/xs_graphite_pure.data")
chiPhysicsTransportXSSet(xs_2,PDT_XSFILE,"CHI_TEST/xs_3_170.data")
chiPhysicsTransportXSSet(xs_3,PDT_XSFILE,"CHI_TEST/xs_air50RH.data")

combo ={{xs_1, 0.5},
        {xs_2, 0.4},
        {xs_3, 0.1}}
aerated_graphite = chiPhysicsTransportXSMakeCombined(combo)


chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              EXISTING,aerated_graphite)
\endcode
 *
 * \return Returns a handle to another cross-section object that contains the
 *         desired combination.
 *
 * \ingroup LuaPhysicsMaterials
 * */
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

  int table_len = lua_rawlen(L,1);

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
  auto new_xs = new chi_physics::TransportCrossSections;

  new_xs->MakeCombined(combinations);

  chi_physics_handler.trnsprt_xs_stack.push_back(new_xs);
  lua_pushnumber(L,chi_physics_handler.trnsprt_xs_stack.size()-1);

  return 1;
}