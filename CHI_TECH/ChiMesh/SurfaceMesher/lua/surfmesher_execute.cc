#include "../../../ChiLua/chi_lua.h"
#include <iostream>
#include "../Predefined/surfmesher_predefined.h"

#include "../../MeshHandler/chi_meshhandler.h"

#include <chi_log.h>
extern ChiLog& chi_log;



//#############################################################################
/** Executes the surface meshing pipeline.

\param ExportLoadBalance bool Optional flag indicating whether to write
                              xy-partition load factors to log. Default=false

\ingroup LuaSurfaceMesher
\author Jan*/
int chiSurfaceMesherExecute(lua_State *L)
{
  int numArgs = lua_gettop(L);

  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();
  chi_log.Log(LOG_ALLVERBOSE_2) << "Executing surface mesher\n";
  cur_hndlr->surface_mesher->Execute();

  bool export_load_balance = false;
  if (numArgs==1) export_load_balance = lua_toboolean(L,1);

  if (export_load_balance)
    cur_hndlr->surface_mesher->PrintLoadBalanceInfo();

  chi_log.Log(LOG_ALLVERBOSE_2)
    << "chiSurfaceMesherExecute: Surface mesher execution completed."
    << std::endl;

  return 0;
}