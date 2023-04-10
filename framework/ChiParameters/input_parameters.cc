#include "input_parameters.h"

#include "parameter_block.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <sstream>
#include <algorithm>
#include <utility>

#define ThrowInputError                                                        \
  throw std::invalid_argument("Input error: " + ObjectType() + "\n" +          \
                              err_stream.str())

#define ExceptionParamNotPresent(param_name)                                   \
  throw std::logic_error(std::string(__PRETTY_FUNCTION__) + ": Parameter \"" + \
                         param_name + "\" not present in list of parameters.")

namespace chi_objects
{

// #################################################################
/**Sets the object type string for more descriptive error messages.*/
void InputParameters::SetObjectType(const std::string& obj_type)
{
  class_name_ = obj_type;
}

// #################################################################
/**Returns the object type string.*/
std::string InputParameters::ObjectType() const { return class_name_; }

// #################################################################
/**Specialization for block type parameters.*/
void InputParameters::AddOptionalParameterBlock(
  const std::string& name,
  ParameterBlock block,
  const std::string& grouping_name)
{
  AddParameter(std::move(block));
  parameter_tags_[name] = InputParameterTag::OPTIONAL;
  parameter_group_names_[name] = grouping_name;
}

// #################################################################
/**Specialization for block type parameters.*/
void InputParameters::AddRequiredParameterBlock(
  const std::string& name, const std::string& grouping_name)
{
  ParameterBlock new_block(name);
  AddParameter(new_block);
  parameter_tags_[name] = InputParameterTag::REQUIRED;
  parameter_group_names_[name] = grouping_name;
}

// #################################################################
/**Returns an error string if the check has not passed. This method
 * first checks whether all the required parameters are supplied. Then
 * it checks that all the parameters supplied actually maps to valid parameters.
 * */
void InputParameters::AssignParameters(const ParameterBlock& params)
{
  std::stringstream err_stream;

  // ================================== Check required parameters
  for (const auto& [param_index, tag] : parameter_tags_)
  {
    if (tag != InputParameterTag::REQUIRED) continue;

    const auto& req_param = GetParam(param_index);
    const auto& req_param_name = req_param.Name();

    if (deprecation_warning_tags_.count(req_param_name) > 0 or
        deprecation_error_tags_.count(req_param_name) > 0 or
        renamed_error_tags_.count(req_param_name) > 0)
      continue;

    if (not params.Has(req_param.Name()))
      err_stream << "Required param \"" << req_param.Name()
                 << "\" not supplied.\n";
  }

  if (not err_stream.str().empty()) ThrowInputError;

  // ================================== Check unused parameters
  {
    for (const auto& param : params.Parameters())
    {
      const auto& param_name = param.Name();
      if (not this->Has(param_name))
        err_stream << "Invalid param \"" << param_name << "\" supplied.\n";
      else if (renamed_error_tags_.count(param_name) > 0)
      {
        err_stream << "Invalid param \"" << param_name << "\" supplied. ";
        err_stream << "The parameter has been renamed. ";
        err_stream << renamed_error_tags_.at(param_name);
      }
    }

    if (not err_stream.str().empty()) ThrowInputError;
  }

  // ================================== Check deprecation warnings
  {
    const auto& dep_warns = deprecation_warning_tags_;
    for (const auto& param : params.Parameters())
    {
      const auto& param_name = param.Name();

      if (this->Has(param_name) and (dep_warns.count(param_name) > 0))
        chi::log.Log0Warning()
          << "Parameter \"" << param_name << "\" has been deprecated "
          << "and will be removed soon.\n"
          << dep_warns.at(param_name);
    }
  }

  // ================================== Check deprecation errors
  {
    const auto& dep_errs = deprecation_error_tags_;
    for (const auto& param : params.Parameters())
    {
      const auto& param_name = param.Name();

      if (this->Has(param_name) and (dep_errs.count(param_name) > 0))
      {
        chi::log.Log0Error()
          << "Parameter \"" << param_name << "\" has been deprecated.\n"
          << dep_errs.at(param_name);
        chi::Exit(EXIT_FAILURE);
      }
    }
  }

  // ================================== Now attempt to assign values
  for (auto& param : params.Parameters())
  {
    auto& input_param = GetParam(param.Name());

    // ====================== Check types match
    if (param.Type() != input_param.Type())
    {
      err_stream << "Invalid parameter type \""
                 << ParameterBlockTypeName(param.Type())
                 << "\" for parameter \"" << param.Name()
                 << "\". Expecting type \""
                 << ParameterBlockTypeName(input_param.Type()) << "\".\n";
      continue;
    } // if type mismatch

    // ====================== Check constraint
    if (constraint_tags_.count(input_param.Name()) != 0)
    {
      const auto& constraint = constraint_tags_.at(input_param.Name());
      if (not constraint->IsAllowable(param.Value()))
      {
        err_stream << constraint->OutOfRangeString(input_param.Name(),
                                                   param.Value());
        err_stream << "\n";
      }
      continue;
    } // if constraint

    input_param = param;
  } // for input params

  if (not err_stream.str().empty()) ThrowInputError;
}

// ##################################################################
/**Marks a parameters as deprecated but will only produce a warning.*/
void InputParameters::MarkParamaterDeprecatedWarning(
  const std::string& param_name, const std::string& deprecation_message /*=""*/)
{
  if (Has(param_name))
    deprecation_warning_tags_[param_name] = deprecation_message;
  else
    ExceptionParamNotPresent(param_name);
}

// ##################################################################
/**Marks a parameters as deprecated and will produce an error if the parameter
 * is specified.*/
void InputParameters::MarkParamaterDeprecatedError(
  const std::string& param_name, const std::string& deprecation_message /*=""*/)
{
  if (Has(param_name))
    deprecation_error_tags_[param_name] = deprecation_message;
  else
    ExceptionParamNotPresent(param_name);
}

// ##################################################################
/**Marks a parameters as renamed and will produce an error if the parameter
 * is specified.*/
void InputParameters::MarkParamaterRenamed(
  const std::string& param_name, const std::string& renaming_description)
{
  if (Has(param_name)) renamed_error_tags_[param_name] = renaming_description;
  else
    ExceptionParamNotPresent(param_name);
}

// ##################################################################
/**Creates a range based constraint for a given parameter.*/
void InputParameters::ConstrainParameterRange(const std::string& param_name,
                                              AllowableRangePtr allowable_range)
{
  if (Has(param_name))
    constraint_tags_[param_name] = std::move(allowable_range);
  else
    ExceptionParamNotPresent(param_name);
}

} // namespace chi_objects