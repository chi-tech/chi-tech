#include "ChiLua/chi_lua.h"
#include<iostream>
#include "chi_runtime.h"

#include "ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h"
#include "ChiPhysics/PhysicsMaterial/material_property_scalarvalue.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"
#include "ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h"


#include "chi_runtime.h"
#include "chi_log.h"
;



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
  const std::string fname = __FUNCTION__;
  const int numArgs = lua_gettop(L);

  if (!((numArgs>=2) && (numArgs<=3)))
  {
    chi::log.Log0Error() << "Incorrect amount of arguments "
                               "in chiPhysicsMaterialAddProperty";
    chi::Exit(EXIT_FAILURE);
  }

  int material_index = lua_tonumber(L,1);
  int property_index = lua_tonumber(L,2);

  const char* provided_name = "";
  if (numArgs == 3)
  {
    provided_name = lua_tostring(L,3);
  }

  //============================================= Get reference to material
  auto cur_material = chi::GetStackItemPtr(chi::material_stack,
                                           material_index, fname);

  //============================================= Process property
  using MatProperty = chi_physics::PropertyType;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCALAR_VALUE
  if (property_index == static_cast<int>(MatProperty::SCALAR_VALUE))
  {
    //Duplicates are allowed

    auto prop = std::make_shared<chi_physics::ScalarValue>();

    prop->property_name = std::string("Property ") +
                          std::to_string(cur_material->properties.size());

    if (numArgs == 3)
      prop->property_name = std::string(provided_name);

    cur_material->properties.push_back(prop);
    chi::log.Log0Verbose1() << "Scalar Value Property added to material"
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
        chi::log.Log0Error()    << "Material " << material_index << " \""
                                   << cur_material->name << "\""
                                   << " already has property "
                                      "TRANSPORT_XSECTIONS"
                                   << std::endl;
        chi::Exit(EXIT_FAILURE);
      }
    }

    auto prop = std::make_shared<chi_physics::TransportCrossSections>();

    prop->property_name = std::string("Property ") +
                          std::to_string(cur_material->properties.size());

    if (numArgs == 3)
      prop->property_name = std::string(provided_name);

    cur_material->properties.push_back(prop);
    chi::log.Log0Verbose1() << "Transport cross-sections added to material"
                                 " at index " << material_index;

    chi::trnsprt_xs_stack.push_back(prop);

    const size_t index = chi::trnsprt_xs_stack.size()-1;

    lua_pushnumber(L,static_cast<lua_Number>(index));
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
        chi::log.Log0Error()    << "Material " << material_index << " \""
                                   << cur_material->name << "\""
                                   << " already has property "
                                      "ISOTROPIC_MG_SOURCE "
                                   << property_index
                                   << std::endl;
        chi::Exit(EXIT_FAILURE);
      }
    }

    auto prop = std::make_shared<chi_physics::IsotropicMultiGrpSource>();

    prop->property_name = std::string("Property ") +
                          std::to_string(cur_material->properties.size());

    if (numArgs == 3)
      prop->property_name = std::string(provided_name);

    cur_material->properties.push_back(prop);
    chi::log.Log0Verbose1() << "Isotropic Multigroup Source added to material"
                                 " at index " << material_index;
  }
  else
  {
    chi::log.Log0Error()
      << "Unsupported property type in call to chiPhysicsMaterialAddProperty.";
    chi::Exit(EXIT_FAILURE);
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
  const std::string fname = __FUNCTION__;
  const int numArgs = lua_gettop(L);

  if (numArgs<3)
  {
    chi::log.Log0Error() << "Incorrect amount of arguments "
                               "in chiPhysicsMaterialSetProperty";
    chi::Exit(EXIT_FAILURE);
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
  auto cur_material = chi::GetStackItemPtr(chi::material_stack,
                                           material_index, fname);

  //============================================= If user supplied name then
  //                                              find property index
  if (!lua_isnumber(L,2))
  {
    for (auto& property : cur_material->properties)
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
      auto prop = std::static_pointer_cast<chi_physics::ScalarValue>(
        cur_material->properties[location_of_prop]);

      //========================== Process operation
      if (operation_index == static_cast<int>(OpType::SINGLE_VALUE))
      {
        double value = lua_tonumber(L,4);
        prop->value = value;
        chi::log.Log0Verbose1() << "Scalar value for material"
                                     " at index " << material_index
                                  << " set to " << value;
      }
      else
      {
        chi::log.Log0Error() << "ERROR: Unsupported operation for "
                                   "SCALAR_VALUE." << std::endl;
        chi::Exit(EXIT_FAILURE);
      }

    }
    else
    {
      chi::log.Log0Error() << "ERROR: Material has no property "
                                 "SCALAR_VALUE." << std::endl;
      chi::Exit(EXIT_FAILURE);
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
      auto prop = std::static_pointer_cast<chi_physics::TransportCrossSections>(
                  cur_material->properties[location_of_prop]);

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

        LuaCheckNilValue("chiPhysicsMaterialSetProperty",L,4);
        int handle = lua_tonumber(L,4);

        std::shared_ptr<chi_physics::TransportCrossSections> xs;
        try {
          xs = chi::GetStackItemPtr(chi::trnsprt_xs_stack, handle, fname);
        }
        catch(const std::out_of_range& o){
          chi::log.LogAllError()
            << "ERROR: Invalid cross-section handle"
            << " in call to chiPhysicsMaterialSetProperty."
            << std::endl;
          chi::Exit(EXIT_FAILURE);
        }
//        auto old_prop = prop;
        prop = xs;

        cur_material->properties[location_of_prop] = prop;

//        delete old_prop; //Still debating if this should be deleted
      }
      else
      {
        chi::log.LogAllError() << "Unsupported operation for "
                                   "TRANSPORT_XSECTIONS." << std::endl;
        chi::Exit(EXIT_FAILURE);
      }

    }
    else
    {
      chi::log.LogAllError() << "Material has no property "
                                 "TRANSPORT_XSECTIONS." << std::endl;
      chi::Exit(EXIT_FAILURE);
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
      auto prop = std::static_pointer_cast<chi_physics::IsotropicMultiGrpSource>(
          cur_material->properties[location_of_prop]);


      if (operation_index == static_cast<int>(OpType::SINGLE_VALUE))
      {
        if (numArgs!=4)
          LuaPostArgAmountError("chiPhysicsMaterialSetProperty",4,numArgs);

        double value = lua_tonumber(L,4);

        prop->source_value_g.resize(1,value);
        chi::log.Log0Verbose1() << "Isotropic Multigroup Source value "
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
          chi::log.LogAllError()
            << "In call to chiPhysicsMaterialSetProperty: "
            << "Material \"" << cur_material->name
            << "\", when setting "
            << "ISOTROPIC_MG_SOURCE using operation FROM_ARRAY, the fourth "
               "argument was detected not to be a lua table.";
          chi::Exit(EXIT_FAILURE);
        }

        const size_t table_len = lua_rawlen(L,4);

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
        chi::log.Log0Verbose1() << "Isotropic Multigroup Source populated "
                                  << " with " << table_len << " values";
      }
      else
      {
        chi::log.LogAllError() << "Unsupported operation for "
                                     "ISOTROPIC_MG_SOURCE." << std::endl;
        chi::Exit(EXIT_FAILURE);
      }
    }
    else
    {
      chi::log.LogAllError() << "Material \"" << cur_material->name
                                << "\" has no property "
                                   "ISOTROPIC_MG_SOURCE." << std::endl;
      chi::Exit(EXIT_FAILURE);
    }
  }
  else
  {
    chi::log.LogAllError() << "Unsupported material property specified in "
                               "call to chiPhysicsMaterialSetProperty."
                               << property_index
                               << std::endl;
    chi::Exit(EXIT_FAILURE);
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
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
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
  auto cur_material = chi::GetStackItemPtr(chi::material_stack,
                                           material_index, fname);

  //============================================= If user supplied name then
  //                                              find property index
  if (!lua_isnumber(L,2))
  {
    for (auto& property : cur_material->properties)
      if (property->property_name == property_index_name)
        property_index = static_cast<int>(property->Type());
  }

  //============================================= Process property
  bool property_polulated = false;
  for (auto& property : cur_material->properties)
  {
    if (static_cast<int>(property->Type()) == property_index)
    {
      property->PushLuaTable(L);
      property_polulated = true;
    }
  }


  if (not property_polulated)
  {
    chi::log.LogAllError() << "Invalid material property specified in "
                                 "call to chiPhysicsMaterialGetProperty."
                              << property_index
                              << std::endl;
    chi::Exit(EXIT_FAILURE);
  }

  return 1;
}


//############################################################
/**Set a cross section to a new value.
\ingroup LuaPhysicsMaterials
 */
int chiPhysicsMaterialModifyTotalCrossSection(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 3)
    LuaPostArgAmountError(__FUNCTION__, 3, num_args);
  LuaCheckNilValue(__FUNCTION__, L, 1);
  LuaCheckNilValue(__FUNCTION__, L, 2);
  LuaCheckNilValue(__FUNCTION__, L, 3);

  int material_index = lua_tonumber(L, 1);
  int group_num = lua_tointeger(L, 2);
  double xs_val = lua_tonumber(L, 3);

  auto cur_material = chi::GetStackItemPtr(chi::material_stack,
                                           material_index, fname);

  using MatProperty = chi_physics::PropertyType;
  using XS = chi_physics::TransportCrossSections;
  bool xs_property_found = false;
  for (auto& property : cur_material->properties)
  {
    auto ptype = MatProperty::TRANSPORT_XSECTIONS;
    if (property->Type() == ptype) {
      auto xs = std::dynamic_pointer_cast<XS>(property);
      xs->sigma_t[group_num] = xs_val;
      chi::log.Log()
        << "sigma_t for group " << group_num
        << " in material " << material_index
        << " changed to " << xs_val;
      xs_property_found = true;
    }
  }
  if (not xs_property_found)
    chi::log.LogAllError()
      << __FUNCTION__ << ": No cross section property for "
                         "specified material.";
  return 0;
}

