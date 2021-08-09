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
private:
  std::string text_name;
public:
  std::vector<chi_mesh::Region*>           regions;
  std::vector<std::shared_ptr<FieldFunction>> field_functions;

public:
  explicit
  Solver(const std::string& in_text_name) : text_name(in_text_name) {}
  virtual ~Solver() {};

  std::string TextName() const {return text_name;}

  //01
  void AddRegion(chi_mesh::Region* region);

  //02
  virtual void Execute() = 0;
};



#endif