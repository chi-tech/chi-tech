#ifndef CHITECH_CHI_OBJECT_H
#define CHITECH_CHI_OBJECT_H

#include "chi_runtime.h"
#include "ChiParameters/input_parameters.h"

class ChiObject
{
private:
  size_t stack_id_ = chi::SIZE_T_INVALID;
  chi_objects::ParameterBlock param_block_used_at_construction_;

public:
  /**Returns the input parameters. For the base ChiObject, there
  * are now parameters loaded.*/
  static chi_objects::InputParameters GetInputParameters();

  /**Default constructor. This will be removed in future.*/
  ChiObject();

  /**Constructor with input parameters.*/
  explicit ChiObject(const chi_objects::InputParameters& params);

  // Setters
  /**Sets the stack id of the object. This allows this
   * object to know its place in the global space.*/
  void SetStackID(size_t stack_id);

  /**Copies the parameter block used to set the input parameters
   * of this object. This is useful to know if an optional parameter
   * was activated/used/set.*/
  void
  SetParamBlockUsedAtConstruction(const chi_objects::ParameterBlock& params);

  // Getters
  /**Returns the stack id of this object. This can be used
   * with input language to connect objects together.*/
  size_t StackID() const;

  /**Sets the parameter block used at construction which allows developers
  * to check for supplied parameters.*/
  const chi_objects::ParameterBlock& ParamBlockUsedAtConstruction() const;

  /**An overridable callback that is called by the ObjectMaker and by default
  * adds the object onto the object stack. This function can be used to
  * place the object on a different stack.*/
  virtual void PushOntoStack(std::shared_ptr<ChiObject>& new_object);

  virtual ~ChiObject() = default;
};

#endif // CHITECH_CHI_OBJECT_H
