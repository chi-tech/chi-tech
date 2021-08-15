#ifndef CHI_BASIC_OPTIONS_H
#define CHI_BASIC_OPTIONS_H

#include "ChiDataTypes/varying.h"

namespace chi_physics
{
//###################################################################
/**Class for option.*/
class BasicOption
{
private:
  std::string     name;
  chi_data_types::Varying         value;

public:
  BasicOption(const std::string& name, const std::string& string_value) :
    name(name), value(string_value) {}

  BasicOption(const std::string& name, const bool& bool_value) :
    name(name), value(bool_value) {}

  BasicOption(const std::string& name, const int64_t& integer_value) :
    name(name), value(integer_value) {}

  BasicOption(const std::string& name, const double& float_value) :
    name(name), value(float_value) {}

  chi_data_types::VaryingDataType Type() const {return value.Type();}

  std::string Name() const {return name;}
  std::string StringValue() const  {return value.StringValue();}
  bool        BoolValue() const    {return value.BoolValue();}
  int64_t     IntegerValue() const {return value.IntegerValue();}
  double      FloatValue() const   {return value.FloatValue();}

  void SetStringValue(const std::string& in_string_value) {value = in_string_value;}
  void SetBoolValue(const bool& in_bool_value)            {value = in_bool_value;}
  void SetIntegerValue(const int64_t& in_integer_value)   {value = in_integer_value;}
  void SetFloatValue(const double& in_float_value)        {value = in_float_value;}

};

//###################################################################
/**Class for basic options*/
class BasicOptions
{
private:
  std::vector<BasicOption> options;

public:
  BasicOptions() = default;

  /**Constructor with initializer list.*/
  BasicOptions(std::initializer_list<BasicOption> in_options) :
    options(in_options)
  {  }

  /**Returns a constant reference to an option that matches the
   * requested name. If no name-match is found the method will throw
   * a std::out_of_range exception.*/
  const BasicOption& operator()(const std::string& option_name) const;

  /**Returns a constant reference to an option at the given
   * index. If the index is out of range then a std::out_of_range
   * exception is thrown. This method can potentially be faster than
   * the string comparison equivalent.*/
  const BasicOption& operator()(size_t index) const;

  /**Returns a non-constant reference to an option that matches the
   * requested name. If no name-match is found the method will throw
   * a std::out_of_range exception.*/
  BasicOption& operator[](const std::string& option_name);

  /**Returns a non-constant reference to an option at the given
   * index. If the index is out of range then a std::out_of_range
   * exception is thrown. This method can potentially be faster than
   * the string comparison equivalent.*/
  BasicOption& operator[](size_t index);

  /**Attempts to find an option that matches the requested name.
   * If one is found then its corresponding index is
   * returned. If it is not found then a std::out_of_range
   * exception is thrown.*/
  size_t GetOptionIndexFromName(const std::string& option_name) const;

};

}//namespace chi_physics
#endif