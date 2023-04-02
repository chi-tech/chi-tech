#ifndef CHITECH_PARAMETER_BLOCK_H
#define CHITECH_PARAMETER_BLOCK_H

#include "varying.h"

#include <memory>
#include <vector>
#include <string>

namespace chi_lua
{
  class TableParserAsParameterBlock;
}

namespace chi_data_types
{

enum class ParameterBlockType
{
  NONE    = -1    /*LUA_TNONE   */,
  Nil     = 0     /*LUA_TNIL    */,
  Boolean = 1     /*LUA_TBOOLEAN*/,
  Number  = 3     /*LUA_TNUMBER */,
  String  = 4     /*LUA_TSTRING */,
  Integer = 5,
  Array   = 98,
  Block  = 99
};

std::string ParameterBlockTypeName(ParameterBlockType type);

class ParameterBlock;
typedef std::unique_ptr<ParameterBlock> ParameterBlockPtr;

class ParameterBlock
{
  friend class chi_lua::TableParserAsParameterBlock;
private:
  ParameterBlockType type_ = ParameterBlockType::NONE;
  std::string keyword_;
  std::unique_ptr<Varying> value_ptr_ = nullptr;
  std::vector<ParameterBlockPtr> parameters_;

public:
  //Helpers
  template<typename T> struct IsBool
  {static constexpr bool value = std::is_same_v<T,bool>;};
  template<typename T> struct IsFloat
  {static constexpr bool value = std::is_floating_point_v<T>;};
  template<typename T> struct IsString
  {static constexpr bool value = std::is_same_v<T,std::string>;};
  template<typename T> struct IsInteger
  {static constexpr bool value = std::is_integral_v<T> and
                                 not std::is_same_v<T,bool>;};

  //Constructors
  /**Constructs an empty parameter block.*/
  explicit ParameterBlock(const std::string& key_str_name = "");

  /**Constructs one of the fundamental types.*/
  template<typename T>
  explicit ParameterBlock(const std::string& key_str_name, T value) :
    keyword_(key_str_name)
  {
    constexpr bool is_supported =
      IsBool<T>::value or IsFloat<T>::value or IsString<T>::value or
      IsInteger<T>::value;

    static_assert(is_supported, "Value type not supported for parameter block");

    if (IsBool<T>::value) type_ = ParameterBlockType::Boolean;
    if (IsFloat<T>::value) type_ = ParameterBlockType::Number;
    if (IsString<T>::value) type_ = ParameterBlockType::String;
    if (IsInteger<T>::value) type_ = ParameterBlockType::Integer;

    value_ptr_ = std::make_unique<Varying>(value);
  }

  //Accessors
  ParameterBlockType Type() const;
  std::string Name() const;
  const Varying& Value() const;
  size_t NumParameters() const;

  //Mutators
protected:
  void ChangeToArray();
public:
  //utilities
protected:
  void AddParameter(ParameterBlockPtr block);
public:
  bool Has(const std::string& param_name) const;

  const ParameterBlock& GetParam(const std::string& param_name) const;
  const ParameterBlock& GetParam(size_t index) const;

public:
  /**Returns the value of the parameter.*/
  template<typename T>
  T GetValue() const
  {
    if (value_ptr_ == nullptr)
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                             ": Value not available for block type " +
                             ParameterBlockTypeName(Type()));
    return Value().GetValue<T>();
  }

  /**Fetches the parameter with the given name and returns it value.*/
  template<typename T>
  T GetParamValue(const std::string& param_name) const
  {
    const auto& param = GetParam(param_name);
    return param.GetValue<T>();
  }

  /**Converts the parameters of an array-type parameter block to a vector of
   * primitive types and returns it.*/
  template<typename T>
  std::vector<T> GetVectorValue() const
  {
    if (Type() != ParameterBlockType::Array)
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                             ": Invalid type requested for parameter of type " +
                             ParameterBlockTypeName(Type()));

    std::vector<T> vec;
    if (parameters_.empty()) return vec;

    //Check the first sub-param is of the right type
    const auto& front_param = parameters_.front();

    //Check that all other parameters are of the required type
    for (const auto& param : parameters_)
      if (param->Type() != front_param->Type())
        throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
          ": Cannot construct vector from block because "
          "the sub_parameters do not all have the correct type.");

    const size_t num_params = parameters_.size();
    for (size_t k=0; k<num_params; ++k)
    {
      const auto& param = GetParam(k);
      vec.push_back(param.GetValue<T>());
    }

    return vec;
  }

  /**Gets a vector of primitive types from an array-type parameter block
   * specified as a parameter of the current block.*/
  template<typename T>
  std::vector<T> GetParamVectorValue(const std::string& param_name) const
  {
    const auto& param = GetParam(param_name);
    return param.GetVectorValue<T>();
  }

  /**Given a reference to a string, recursively travels the parameter
   * tree and print values into the reference string.*/
  void RecursiveDumpToString(std::string& outstr,
                             const std::string& offset="") const;
};

}//namespace chi_lua

#endif //CHITECH_PARAMETER_BLOCK_H
