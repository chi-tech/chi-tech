#include "../../ChiLua/chi_lua.h"
#include<iostream>
#include "../chi_physics.h"
#include "../PhysicsMaterial/chi_physicsmaterial.h"

extern ChiPhysics chi_physics_handler;

#include <chi_log.h>

extern ChiLog chi_log;

/** \defgroup LuaPhysicsMaterials Materials
 * \ingroup LuaPhysics*/

//#############################################################################
/** Adds a material to the problem. Materials are added to the global
 * physics handler and is therefore accessible across all meshes and solvers.
 *
\param Name char (Optional) Material name.

\return MaterialHandle int Handle to the created material.


##_

###Example\n
Simple example of adding a material
\code
materials = {}
materials[0] = chiPhysicsAddMaterial("Test Material");
\endcode


\ingroup LuaPhysicsMaterials
\author Jan*/
int chiPhysicsAddMaterial(lua_State *L)
{
  int numArgs = lua_gettop(L);

  chi_physics::Material* new_material = new chi_physics::Material;
  if (numArgs==1)
  {
    const char* temp = lua_tostring(L, 1);
    new_material->name = std::string(temp);
  }

  chi_physics_handler.material_stack.push_back(new_material);

  int index = chi_physics_handler.material_stack.size()-1;
  lua_pushnumber(L,index);

  chi_log.Log(LOG_0VERBOSE_1) << "New material added at index " << index
                            << " with name \"" << new_material->name << "\"";

  return 1;
}
