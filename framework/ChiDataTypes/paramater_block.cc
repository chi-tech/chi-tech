#include "parameter_block.h"

#include <cmath>
#include <algorithm>

namespace chi_data_types
{
//#################################################################
std::string ParameterBlockTypeName(ParameterBlockType type)
{
  switch (type)
  {
    case ParameterBlockType::NONE:    return "NONE";
    case ParameterBlockType::Nil:     return "NIL";
    case ParameterBlockType::Boolean: return "BOOLEAN";
    case ParameterBlockType::Number:  return "NUMBER";
    case ParameterBlockType::String:  return "STRING";
    case ParameterBlockType::Integer: return "INTEGER";
    case ParameterBlockType::Array:   return "ARRAY";
    case ParameterBlockType::Block:   return "BLOCK";
    default:
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
      ": No name associated with type");
  }
}

//#################################################################
ParameterBlock::ParameterBlock(const std::string& key_str_name) :
  type_(ParameterBlockType::Block),
  keyword_(key_str_name)
{}

//Accessors
ParameterBlockType ParameterBlock::Type() const {return type_;}
std::string ParameterBlock::Name() const {return keyword_;}

//#################################################################
const Varying& ParameterBlock::Value() const
{
  switch (this->Type())
  {
    case ParameterBlockType::Boolean:
    case ParameterBlockType::Number:
    case ParameterBlockType::String:
    case ParameterBlockType::Integer:
    {
      if (value_ptr_ == nullptr)
        throw std::runtime_error(std::string(__PRETTY_FUNCTION__) +
          ": Uninitialized Varying value for block " + this->Name());
      return *value_ptr_;
    }
    default:
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
      ": Called for block of type " + ParameterBlockTypeName(this->Type()) +
      " which has no value.");
  }
}

//#################################################################
/**Returns the number of parameters in a block. This is normally
 * only useful for the ARRAY type.*/
size_t ParameterBlock::NumParameters() const
{
  return parameters_.size();
}

//Mutators
//#################################################################
/**Changes the block type to array, making it accessible via integer
 * keys.*/
void ParameterBlock::ChangeToArray()
{
  type_ = ParameterBlockType::Array;
}

//#################################################################
void ParameterBlock::AddParameter(ParameterBlockPtr block)
{
  for (const auto& param : parameters_)
    if (param->Name() == block->Name())
      throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
      ": Attempting to add duplicate parameter " + param->Name() + " to "
      "block " + this->Name());
  parameters_.push_back(std::move(block));

  SortParameters();
}

//#################################################################
void ParameterBlock::SortParameters()
{
  struct AlphabeticFunctor
  {
    bool operator()(const ParameterBlockPtr& paramA,
                    const ParameterBlockPtr& paramB)
    {
      return paramA->Name() < paramB->Name();
    }
  };

  std::sort(parameters_.begin(), parameters_.end(), AlphabeticFunctor());
}

//#################################################################
bool ParameterBlock::Has(const std::string &param_name) const
{
  for (const auto& param : parameters_)
    if (param->Name() == param_name)
      return true;

  return false;
}

//#################################################################
/**Gets a parameter by name.*/
const ParameterBlock& ParameterBlock::
  GetParam(const std::string &param_name) const
{
  for (const auto& param : parameters_)
    if (param->Name() == param_name)
      return *param;

  throw std::out_of_range(std::string(__PRETTY_FUNCTION__) +
  ": Parameter " + param_name + " not present in block");
}

//#################################################################
/**Gets a parameter by index.*/
const ParameterBlock& ParameterBlock::GetParam(size_t index) const
{
  try
  {
    return *parameters_.at(index);
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range(std::string(__PRETTY_FUNCTION__) +
      ": Parameter with index " + std::to_string(index) +
      " not present in block");
  }
}


//#################################################################
/**Print the block tree structure into a designated string.*/
void ParameterBlock::RecursiveDumpToString(std::string& outstr,
                                           const std::string& offset) const
{
  outstr += offset + this->Name() + " = \n";
  outstr += offset + "{\n";

  for (const auto& param : parameters_)
  {

    switch (param->Type())
    {
      case ParameterBlockType::Boolean:
      {
        outstr += offset + "  " + param->Name() + " = ";
        const bool value = param->Value().BoolValue();
        outstr += std::string( value ? "true" : "false") + ",\n";
        break;
      }
      case ParameterBlockType::Number:
      {
        outstr += offset + "  " + param->Name() + " = ";
        const double value = param->Value().FloatValue();
        outstr += std::to_string(value) + ",\n";
        break;
      }
      case ParameterBlockType::String:
      {
        outstr += offset + "  " + param->Name() + " = ";
        const auto& value = param->Value().StringValue();
        outstr += "\"" + value + "\",\n";
        break;
      }
      case ParameterBlockType::Integer:
      {
        outstr += offset + "  " + param->Name() + " = ";
        const int64_t value = param->Value().IntegerValue();
        outstr += std::to_string(value) + ",\n";
        break;
      }
      case ParameterBlockType::Array:
      case ParameterBlockType::Block:
      {
        param->RecursiveDumpToString(outstr, offset + "  ");
        break;
      }
      default: break;
    }
  }//for parameter

  outstr += offset + "}\n";
}

}//namespace chi_data_types