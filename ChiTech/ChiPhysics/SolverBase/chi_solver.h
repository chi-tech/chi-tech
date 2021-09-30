#ifndef CHI_PHYSICS_SOLVER_H
#define CHI_PHYSICS_SOLVER_H
#include<iostream>
#include "../chi_physics_namespace.h"
#include "../../ChiMesh/Region/chi_region.h"

#include "ChiPhysics/BasicOptions/basic_options.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"


/**\defgroup LuaSolver Solvers
 * \ingroup LuaPhysics*/

//######################################################### Solver parent class
class chi_physics::Solver
{
private:
  std::string text_name;
public:
  BasicOptions basic_options;
  std::vector<chi_mesh::Region*> regions;
  std::vector<std::shared_ptr<FieldFunction>> field_functions;

public:
  explicit
  Solver(const std::string& in_text_name) : text_name(in_text_name) {}
  Solver(const std::string& in_text_name,
         std::initializer_list<BasicOption> in_options) :
         text_name(in_text_name),
         basic_options(in_options) {}
  virtual ~Solver() {};

  std::string TextName() const {return text_name;}

  void AddRegion(chi_mesh::Region* region);
  virtual void Initialize() = 0;
  virtual void Execute() = 0;
};



#endif