#ifndef _chi_solver_h
#define _chi_solver_h
#include<iostream>
#include "../chi_physics_namespace.h"
#include "../../ChiMesh/CHI_REGION/chi_region.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"


/**\defgroup LuaSolver Solvers
 * \ingroup LuaPhysics*/

//######################################################### Solver parent class
class chi_physics::Solver
{
public:
  std::vector<chi_mesh::Region*>           regions;
  std::vector<chi_physics::FieldFunction*> field_functions;

public:
  //01
  void AddRegion(chi_mesh::Region* region);

  //02
  virtual void Execute();
};



#endif