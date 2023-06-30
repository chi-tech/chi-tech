#ifndef CHI_PHYSICS_SOLVER_H
#define CHI_PHYSICS_SOLVER_H

#include "ChiObject.h"
#include "physics/chi_physics_namespace.h"

#include "physics/BasicOptions/basic_options.h"
#include"parameters/parameter_block.h"

#include <iostream>
#include <utility>

namespace chi_physics
{
  class FieldFunctionGridBased;
}


/**\defgroup LuaSolver Solvers
 * \ingroup LuaPhysics*/

//######################################################### Solver parent class
class chi_physics::Solver : public ChiObject
{
private:
  const std::string text_name_;
protected:
  BasicOptions basic_options_;
  std::vector<std::shared_ptr<FieldFunctionGridBased>> field_functions_;

public:
  static chi::InputParameters GetInputParameters();
  explicit
  Solver(std::string  in_text_name) : text_name_(std::move(in_text_name)) {}
  Solver(std::string  in_text_name,
         std::initializer_list<BasicOption> in_options) :
    text_name_(std::move(in_text_name)),
    basic_options_(in_options) {}
  explicit Solver(const chi::InputParameters& params);
  virtual ~Solver() = default;

  BasicOptions& GetBasicOptions() {return basic_options_;}
  const BasicOptions& GetBasicOptions() const {return basic_options_;}

  std::vector<std::shared_ptr<FieldFunctionGridBased>>&
  GetFieldFunctions() {return field_functions_;}

  const std::vector<std::shared_ptr<FieldFunctionGridBased>>&
  GetFieldFunctions() const {return field_functions_;}

  std::string TextName() const {return text_name_;}

  virtual void Initialize();
  virtual void Execute();
  virtual void Step();
  virtual void Advance();
};



#endif