#include"../../../CHI_LUA/chi_lua.h"
#include<iostream>
#include "../chi_linemesh.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"

/** \defgroup LuaLineMesh Line Meshes
 * \ingroup LuaMesh
*/

#include <chi_log.h>

extern CHI_LOG chi_log;

//#############################################################################
/** Creates a new line mesh from a loop.
 *
\param LoopCollectionHandle int Handle to the Loop collection.
\param LoopHandle int Handle to the loop inside the collection.

\return Handle int Handle to the created line mesh.
\ingroup LuaLineMesh
\author Jan*/
int chiLineMeshCreateFromLoop(lua_State *L)
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

  //============================================= Check valid loop
  if (cur_loop->edges.size()==0)
  {
    std::cerr << "ERROR: Empty loop accessed." << std::endl;
    exit(EXIT_FAILURE);
  }

  //============================================= Get vertices
  chi_mesh::Vertex vi,vf;

  vi=cur_loop->edges.front().vertices[0];
  vf=cur_loop->edges.back().vertices[1];

  //printf("vi %+.4f %+.4f %+.4f\n",vi.x,vi.y,vi.z);
  //printf("vf %+.4f %+.4f %+.4f\n",vf.x,vf.y,vf.z);

  //============================================= Create LineMesh
  chi_mesh::LineMesh* new_line = new chi_mesh::LineMesh;
  new_line->vertices.push_back(vi);
  new_line->vertices.push_back(vf);

  //============================================= Add to handler
  cur_hndlr->linemesh_stack.push_back(new_line);

  int index = cur_hndlr->linemesh_stack.size()-1;
  lua_pushnumber(L,index);

  chi_log.Log(LOG_ALLVERBOSE_2)
  << "chiLineMeshCreateFromLoop: Created line mesh " << index << std::endl;

  return 1;
}