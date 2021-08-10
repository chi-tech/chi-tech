#ifndef CHI_BASIC_OPTIONS_H
#define CHI_BASIC_OPTIONS_H

namespace chi_physics
{

enum class BasicOptionType : int
{
  STRING  = 1,
  BOOL    = 2,
  INTEGER = 3,
  FLOAT   = 4
};

//###################################################################
/**Class for option.*/
class BasicOption
{
private:
  BasicOptionType type;
  std::string     name;

  std::string     string_value;
  bool            bool_value = false;
  int64_t         integer_value = 0;
  double          float_value = 0.0;

public:
  BasicOption(const std::string& name, const std::string& string_value) :
    type(BasicOptionType::STRING), name(name), string_value(string_value) {}

  BasicOption(const std::string& name, const bool& bool_value) :
    type(BasicOptionType::BOOL), name(name), bool_value(bool_value) {}

  BasicOption(const std::string& name, const int64_t& integer_value) :
    type(BasicOptionType::INTEGER), name(name), integer_value(integer_value) {}

  BasicOption(const std::string& name, const double& float_value) :
    type(BasicOptionType::FLOAT), name(name), float_value(float_value) {}

  BasicOptionType Type() const {return type;}

  std::string Name() const {return name;}
  std::string StringValue() const {return string_value;}
  bool        BoolValue() const {return bool_value;}
  int64_t     IntegerValue() const {return integer_value;}
  double      FloatValue() const {return float_value;}

  void SetStringValue(const std::string& in_string_value) {string_value = in_string_value;}
  void SetBoolValue(const bool& in_bool_value)            {bool_value = in_bool_value;}
  void SetIntegerValue(const int64_t& in_integer_value)   {integer_value = in_integer_value;}
  void SetFloatValue(const double& in_float_value)        {float_value = in_float_value;}

};

//###################################################################
/**Class for basic options*/
class BasicOptions
{
private:
  std::vector<BasicOption> options;

public:
  BasicOptions() = default;

  BasicOptions(std::initializer_list<BasicOption> in_options)
  {
    for (auto& option : in_options)
      options.push_back(option);
  }

  const BasicOption& operator()(const std::string& option_name)
  {
    for (const auto& option : options)
      if (option.Name() == option_name)
        return option;

    throw std::out_of_range("Basic option " + option_name +
                            " does not appear to exist.");
  }

  const BasicOption& operator()(size_t index)
  {
    if (index < options.size())
      return options[index];

    throw std::out_of_range("Basic option with index " + std::to_string(index) +
    " does not appear to exist.");
  }

  BasicOption& operator[](const std::string& option_name)
  {
    for (auto& option : options)
      if (option.Name() == option_name)
        return option;

    throw std::out_of_range("Basic option " + option_name +
    " does not appear to exist.");
  }

  BasicOption& operator[](size_t index)
  {
    if (index < options.size())
      return options[index];

    throw std::out_of_range("Basic option with index " + std::to_string(index) +
    " does not appear to exist.");
  }

};

}//namespace chi_physics
#endif