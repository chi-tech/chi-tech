#ifndef CHI_PHYSICS_SOLVER_H
#define CHI_PHYSICS_SOLVER_H

#include "ChiObject.h"
#include "physics/chi_physics_namespace.h"

#include "physics/BasicOptions/basic_options.h"
#include "parameters/parameter_block.h"

#include <iostream>
#include <utility>

namespace chi_physics
{
class FieldFunctionGridBased;
class TimeStepController;
}

/**\defgroup LuaSolver Solvers
 * \ingroup LuaPhysics*/

// ######################################################### Solver parent class
class chi_physics::Solver : public ChiObject
{
public:
  static chi::InputParameters GetInputParameters();
  explicit Solver(std::string in_text_name)
    : text_name_(std::move(in_text_name))
  {
  }
  Solver(std::string in_text_name,
         std::initializer_list<BasicOption> in_options)
    : text_name_(std::move(in_text_name)), basic_options_(in_options)
  {
  }
  explicit Solver(const chi::InputParameters& params);
  virtual ~Solver() = default;

  std::string TextName() const;

  BasicOptions& GetBasicOptions();
  const BasicOptions& GetBasicOptions() const;

  std::vector<std::shared_ptr<FieldFunctionGridBased>>& GetFieldFunctions();

  const std::vector<std::shared_ptr<FieldFunctionGridBased>>&
  GetFieldFunctions() const;

  virtual double TimeStepSize() const;
  virtual double Time() const;
  virtual double EndTime() const;
  virtual int MaxTimeSteps() const;
  virtual size_t TimeStepIndex() const;

  virtual void Initialize();
  virtual void Execute();
  virtual void Step();
  virtual void Advance();

  /**Generalized query for information supporting varying returns.*/
  virtual chi::ParameterBlock GetInfo(const chi::ParameterBlock& params) const;
  /**PreCheck call to GetInfo.*/
  chi::ParameterBlock
  GetInfoWithPreCheck(const chi::ParameterBlock& params) const;

private:
  const std::string text_name_;

protected:
  BasicOptions basic_options_;
  std::vector<std::shared_ptr<FieldFunctionGridBased>> field_functions_;

  std::shared_ptr<TimeStepController> time_step_controller_ = nullptr;

  double dt_ = 0.01;
  double time_ = 0.0;
  double end_time_ = 1.0;
  int max_time_steps_ = -1;

  size_t t_index_ = 0;
};

#endif