#ifndef CHITECH_INPUT_PARAMETERS_H
#define CHITECH_INPUT_PARAMETERS_H

#include "parameter_block.h"

#include <map>

#define MakeInpParamsForObj(x, y)                                              \
  chi_objects::InputParameters::MakeForObject<x>(y)

namespace chi_objects
{

enum class InputParameterTag
{
  NONE = 0,
  OPTIONAL = 1,
  REQUIRED = 2
};

/**Class for storing input parameters.*/
class InputParameters : public ParameterBlock
{
private:
  std::string class_name_;
  std::map<std::string, InputParameterTag> parameter_tags_;
  std::map<std::string, std::string> parameter_group_names_;

  std::map<std::string, std::string> deprecation_warning_tags_;
  std::map<std::string, std::string> deprecation_error_tags_;
  std::map<std::string, std::string> renamed_error_tags_;

public:
  InputParameters() = default;
  template <typename T>
  static InputParameters MakeForObject(const ParameterBlock& params)
  {
    auto input_param = T::GetInputParameters();

    input_param.AssignParameters(params);
    return input_param;
  }

public:
  void SetObjectType(const std::string& obj_type);
  std::string ObjectType() const;

private:
  using ParameterBlock::AddParameter;

public:
  template <typename T>
  void AddOptionalParameter(const std::string& name,
                            T value,
                            const std::string& grouping_name = "")
  {
    AddParameter(name, value);
    parameter_tags_[name] = InputParameterTag::OPTIONAL;
    parameter_group_names_[name] = grouping_name;
  }
  void AddOptionalParameterBlock(const std::string& name,
                                 ParameterBlock block,
                                 const std::string& grouping_name = "");

  template <typename T>
  void AddRequiredParameter(const std::string& name,
                            const std::string& grouping_name = "")
  {
    AddParameter(name, chi_data_types::Varying::DefaultValue<T>());
    parameter_tags_[name] = InputParameterTag::REQUIRED;
    parameter_group_names_[name] = grouping_name;
  }

  void AddRequiredParameterBlock(const std::string& name,
                                 const std::string& grouping_name = "");

  template <typename T>
  void ChangeExistingParamToOptional(const std::string& name,
                                     T value,
                                     const std::string& grouping_name = "")
  {
    auto& param = GetParam(name);
    param = ParameterBlock(name, value);
    parameter_tags_[name] = InputParameterTag::OPTIONAL;
    parameter_group_names_[name] = grouping_name;
  }

  template <typename T>
  void ChangeExistingParamToRequired(const std::string& name,
                                     const std::string& grouping_name = "")
  {
    auto& param = GetParam(name);
    param = ParameterBlock(name, chi_data_types::Varying::DefaultValue<T>());
    parameter_tags_[name] = InputParameterTag::REQUIRED;
    parameter_group_names_[name] = grouping_name;
  }

public:
  void AssignParameters(const ParameterBlock& params);

  void MarkParamaterDeprecatedWarning(const std::string& param_name,
                                      const std::string& deprecation_message="");
  void MarkParamaterDeprecatedError(const std::string& param_name,
                                    const std::string& deprecation_message="");
  void MarkParamaterRenamed(const std::string& param_name,
                            const std::string& renaming_description);
};

} // namespace chi_objects

#endif // CHITECH_INPUT_PARAMETERS_H
