#include "chi_mesh.h"
#include<iostream>

#include "MeshHandler/chi_meshhandler.h"

extern std::vector<chi_mesh::MeshHandler*>  chi_meshhandler_stack;
extern int                                  chi_current_mesh_handler;

//###################################################################
/**Obtains a pointer to the current mesh handler from the global stack.
\author Jan*/
chi_mesh::MeshHandler* chi_mesh::GetCurrentHandler()
{
  chi_mesh::MeshHandler* cur_handler;

  try{
    cur_handler = chi_meshhandler_stack.at(chi_current_mesh_handler);
  }
  catch(const std::out_of_range& err){
    std::cerr << "ERROR: Invalid index to mesh handler.";
    std::cerr << " Call chiMeshHandlerCreate() to create handler scope.\n";
    exit(EXIT_FAILURE);
  }

  return cur_handler;
}

//###################################################################
/**Adds a new mesh handler to the stack, sets it as the current handler
 * and returns a handle to it.*/
size_t chi_mesh::PushNewHandler()
{
  auto new_handler = new chi_mesh::MeshHandler;

  chi_meshhandler_stack.push_back(new_handler);

  int index = chi_meshhandler_stack.size()-1;
  chi_current_mesh_handler = index;

  return index;
}

//###################################################################
/**Adds a new mesh handler to the stack, sets it as the current handler
 * and returns a pointer to it.*/
chi_mesh::MeshHandler* chi_mesh::GetNewHandler()
{
  auto new_handler = new chi_mesh::MeshHandler;

  chi_meshhandler_stack.push_back(new_handler);

  int index = chi_meshhandler_stack.size()-1;
  chi_current_mesh_handler = index;

  return new_handler;
}


