
#include "ChiConsole/chi_console.h"

#include "ChiMath/lua/chi_math_lua.h"
#include "ChiMesh/lua/chi_mesh_lua.h"
#include "ChiMPI/lua/chi_mpi_lua.h"
#include "ChiLog/lua/chi_log_lua.h"
#include "ChiPhysics/lua/physics_lua_utils.h"
#include "LuaTest/lua_test.h"

#include "lua/chi_modules_lua.h"

//############################################################################# Default constructor
/** Default constructor for the console*/
chi_objects::ChiConsole::ChiConsole() noexcept
{
	//========================================== Initializing console
	consoleState = luaL_newstate();
  auto& L = this->consoleState;

	luaL_openlibs(L);

	//========================================== Registering functions
//	#include"../ChiLua/chi_lua_register.h"

  chi_math::lua_utils::RegisterLuaEntities(L);
  chi_mesh::lua_utils::RegisterLuaEntities(L);
  chi_mpi_utils::lua_utils::RegisterLuaEntities(L);
  chi_log_utils::lua_utils::RegisterLuaEntities(L);
  chi_physics::lua_utils::RegisterLuaEntities(L);
  chi_lua_test::lua_utils::RegisterLuaEntities(L);

  chi_modules::lua_utils::RegisterLuaEntities(L);
}


