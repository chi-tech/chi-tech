#include "../../../CHI_LUA/chi_lua.h"
#include <iostream>
#include "../Predefined2D/volmesher_predefined2d.h"

#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include <chi_log.h>
#include <CHI_TIMER/chi_timer.h>

extern CHI_LOG chi_log;
extern CHI_TIMER chi_program_timer;

#include "../../../CHI_CONSOLE/chi_console.h"
extern CHI_CONSOLE       chi_console;

#include <iomanip>


//#############################################################################
/** Executes the volume meshing pipeline.

\ingroup LuaVolumeMesher
\author Jan*/
int chiVolumeMesherExecute(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //Get memory before
  CST_MEMORY mem_before = chi_console.GetMemoryUsage();

  cur_hndlr->volume_mesher->Execute();

  //Get memory usage
  CST_MEMORY mem_after = chi_console.GetMemoryUsage();

  std::stringstream mem_string;
  mem_string
  << " Memory used = " << std::setprecision(3)
  << mem_after.memory_mbytes - mem_before.memory_mbytes
  << " MB\n"
  << "Total process memory used after meshing "
  << mem_after.memory_mbytes
  << " MB";

  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << "chiVolumeMesherExecute: Volume meshing completed."
    << mem_string.str()
    << std::endl;

  return 0;
}