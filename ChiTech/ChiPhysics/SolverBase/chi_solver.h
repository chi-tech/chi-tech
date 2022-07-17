#ifndef CHI_PHYSICS_SOLVER_H
#define CHI_PHYSICS_SOLVER_H
#include<iostream>
#include <utility>
#include "ChiPhysics/chi_physics_namespace.h"

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
  std::vector<std::shared_ptr<FieldFunction>> field_functions;

public:
  explicit
  Solver(std::string  in_text_name) : text_name(std::move(in_text_name)) {}
  Solver(std::string  in_text_name,
         std::initializer_list<BasicOption> in_options) :
         text_name(std::move(in_text_name)),
         basic_options(in_options) {}
  virtual ~Solver() = default;

  std::string TextName() const {return text_name;}

  virtual void Initialize() = 0;
  virtual void Execute() = 0;
};



#endif