#include "../../CHI_LUA/chi_lua.h"
#include <iostream>
#include <sstream>



#include "../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../chi_mesh_edgeloops.h"

#include <chi_log.h>

extern ChiLog chi_log;

//#############################################################################
/** Splits an edge loop into edges if they differ by a certain angle.

\param LoopCollectionHandle int Handle to the Loop collection.
\param LoopHandle int Handle to the loop inside the collection.
\param Angle double (Optional) Value of the angle by which to split. Default 1 deg.

\return Handle int. Handle to the newly created LoopCollection.
\ingroup LuaMesh
\author Jan*/
int chiEdgeLoopSplitByAngle(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Get the loop collection
  int loop_coll_index = lua_tonumber(L,1);

  chi_mesh::EdgeLoopCollection* cur_coll;
  try{
    cur_coll = cur_hndlr->edge_loop_collections.at(loop_coll_index);
  }
  catch(std::out_of_range o)
  {
    std::cerr << "ERROR: Invalid loop collection handle." << std::endl;
    exit(EXIT_FAILURE);
  }

  //============================================= Get the loop
  int loop_index = lua_tonumber(L,2);

  chi_mesh::EdgeLoop* cur_loop;
  try{
    cur_loop = cur_coll->at(loop_index);
  }
  catch(std::out_of_range o)
  {
    std::cerr << "ERROR: Invalid loop handle." << std::endl;
    exit(EXIT_FAILURE);
  }

  //============================================= Get the angle
  bool angle_specified = false;
  double angle = 1.0;
  if (lua_gettop(L)==3)
  {
    angle = lua_tonumber(L,3);
    angle_specified = true;
  }

  //============================================= Call split function
  chi_mesh::EdgeLoopCollection* new_coll = new chi_mesh::EdgeLoopCollection;
  if (angle_specified)
  {
    new_coll = chi_mesh::SplitEdgeLoopByAngle(cur_loop,angle);
  }
  else
  {
    new_coll = chi_mesh::SplitEdgeLoopByAngle(cur_loop);
  }

  cur_hndlr->edge_loop_collections.push_back(new_coll);

  int index = cur_hndlr->edge_loop_collections.size()-1;
  lua_pushnumber(L,index);
  lua_pushnumber(L,new_coll->size());

  std::stringstream outtext;
  outtext << "chiEdgeLoopSplitByAngle: Number of edge loops after split = ";
  outtext << new_coll->size();
  outtext << std::endl;

  chi_log.Log(LOG_ALLVERBOSE_2) << outtext.str();

  return 2;
}