#include "chi_mesh.h"
#include<iostream>

#include "MeshHandler/chi_meshhandler.h"
#include "chi_runtime.h"


//###################################################################
/**Obtains a reference to the current mesh handler from the global stack.
 *
 * If the stack is empty this routine will through `std::logic_error`.
\author Jan*/
chi_mesh::MeshHandler& chi_mesh::GetCurrentHandler()
{
  if (chi::meshhandler_stack.empty())
    throw std::logic_error("chi_mesh::GetCurrentHandler: No handlers on stack");

  return chi::GetStackItem<chi_mesh::MeshHandler>(chi::meshhandler_stack,
                                                  chi::current_mesh_handler);
}

//###################################################################
/**Adds a new mesh handler to the stack, sets it as the current handler
 * and returns a handle to it.*/
size_t chi_mesh::PushNewHandlerAndGetIndex()
{
  chi::meshhandler_stack.push_back(std::make_shared<chi_mesh::MeshHandler>());

  int index = (int)chi::meshhandler_stack.size() - 1;
  chi::current_mesh_handler = index;

  return index;
}


