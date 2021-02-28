#ifndef CHI_PHYSICS_SOLVER_H
#define CHI_PHYSICS_SOLVER_H
#include<iostream>
#include "../chi_physics_namespace.h"
#include "../../ChiMesh/Region/chi_region.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"


/**\defgroup LuaSolver Solvers
 * \ingroup LuaPhysics*/

//######################################################### Solver parent class
class chi_physics::Solver
{
public:
  std::vector<chi_mesh::Region*>           regions;
  std::vector<std::shared_ptr<FieldFunction>> field_functions;

public:
  virtual ~Solver() {};
  //01
  void AddRegion(chi_mesh::Region* region);

  //02
  virtual void Execute();
};



#endif