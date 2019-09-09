#include "chi_mesh.h"
#include<iostream>

#include "CHI_MESHHANDLER/chi_meshhandler.h"

extern std::vector<chi_mesh::MeshHandler*>  chi_meshhandler_stack;
extern int                                  chi_current_mesh_handler;

/**Obtains a pointer to the current mesh handler from the global stack.
\author Jan*/
chi_mesh::MeshHandler* chi_mesh::GetCurrentHandler()
{
  chi_mesh::MeshHandler* cur_handler;

  try{
    cur_handler = chi_meshhandler_stack.at(chi_current_mesh_handler);
  }
  catch(std::out_of_range err){
    std::cerr << "ERROR: Invalid index to mesh handler.";
    std::cerr << " Call chiMeshHandlerCreate() to create handler scope.\n";
    exit(EXIT_FAILURE);
  }

  return cur_handler;
}