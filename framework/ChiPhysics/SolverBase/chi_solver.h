#ifndef CHI_PHYSICS_SOLVER_H
#define CHI_PHYSICS_SOLVER_H
#include <iostream>
#include <utility>
#include "ChiPhysics/chi_physics_namespace.h"

#include "ChiPhysics/BasicOptions/basic_options.h"
#include "ChiParameters/parameter_block.h"

namespace chi_physics
{
  class FieldFunction;
}


/**\defgroup LuaSolver Solvers
 * \ingroup LuaPhysics*/

//######################################################### Solver parent class
class chi_physics::Solver
{
private:
  std::string text_name_;
protected:
  BasicOptions basic_options_;
  std::vector<std::shared_ptr<FieldFunction>> field_functions_;

public:
  explicit
  Solver(std::string  in_text_name) : text_name_(std::move(in_text_name)) {}
  Solver(std::string  in_text_name,
         std::initializer_list<BasicOption> in_options) :
    text_name_(std::move(in_text_name)),
    basic_options_(in_options) {}
  virtual ~Solver() = default;

  BasicOptions& GetBasicOptions() {return basic_options_;}
  const BasicOptions& GetBasicOptions() const {return basic_options_;}

  std::vector<std::shared_ptr<FieldFunction>>&
  GetFieldFunctions() {return field_functions_;}

  const std::vector<std::shared_ptr<FieldFunction>>&
  GetFieldFunctions() const {return field_functions_;}

  std::string TextName() const {return text_name_;}

  virtual void Initialize();
  virtual void Execute();
  virtual void Step();
};



#endif