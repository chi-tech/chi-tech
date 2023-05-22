#include "varying.h"

#include <algorithm>
#include <sstream>
#include <iomanip>

/**Provides a string-name for an enumerated VaryingDataType.*/
std::string
chi_data_types::VaryingDataTypeStringName(chi_data_types::VaryingDataType type)
{
  switch (type)
  {
    case VaryingDataType::VOID:
      return "VOID";
    case VaryingDataType::ARBITRARY_BYTES:
      return "ARBITRARY_BYTES";
    case VaryingDataType::STRING:
      return "STRING";
    case VaryingDataType::BOOL:
      return "BOOL";
    case VaryingDataType::INTEGER:
      return "INTEGER";
    case VaryingDataType::FLOAT:
      return "FLOAT";
    default:
      return "UNKNOWN";
  }
}

// ###################################################################
/**PopulateRaw template specialization for std::string.*/
template <>
void chi_data_types::Varying::PopulateRaw(const std::string& value)
{
  raw_data_.resize(value.size());
  std::transform(value.begin(),
                 value.end(),
                 raw_data_.begin(),
                 [](char c) { return std::byte(c); });
  if (value.back() != '\0') raw_data_.push_back(std::byte('\0'));

  data_initialized_ = true;
}

/**Checks if two VaryingDataType values match.
 * Type A is matched against type B.*/
void chi_data_types::Varying::CheckTypeMatch(
  const VaryingDataType type_A, const VaryingDataType type_B_required) const
{
  if (type_A != type_B_required)
    throw std::logic_error("Varying data type " + TypeName() +
                           " does not "
                           "correspond to the required type, " +
                           VaryingDataTypeStringName(type_B_required));
}

/**Checks whether the data has been initialized.*/
void chi_data_types::Varying::CheckDataInitialized() const
{
  if (not data_initialized_)
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                           "Varying data type used uninitialized.");
}

// ###################################################################
// Constructors
/**Constructor for an arbitrary sequence of bytes value.*/
chi_data_types::Varying::Varying(const std::vector<std::byte>& value)
  : type_(VaryingDataType::ARBITRARY_BYTES)
{
  raw_data_ = value;
  data_initialized_ = true;
}

/**Constructor for a string value.*/
chi_data_types::Varying::Varying(const std::string& value)
  : type_(VaryingDataType::STRING)
{
  PopulateRaw<std::string>(value);
}

/**Copy constructor.*/
chi_data_types::Varying::Varying(const Varying& other)
{
  raw_data_ = other.raw_data_;
  data_initialized_ = other.data_initialized_;
  type_ = other.type_;
}

/**Move constructor.*/
chi_data_types::Varying::Varying(Varying&& other) noexcept
{
  std::swap(raw_data_, other.raw_data_);
  std::swap(data_initialized_, other.data_initialized_);
  std::swap(type_, other.type_);
}

/**Assignment operator. i.e., type_A = type_B*/
chi_data_types::Varying&
chi_data_types::Varying::operator=(const Varying& other)
{
  if (this != &other)
  {
    raw_data_ = other.raw_data_;
    data_initialized_ = other.data_initialized_;
    type_ = other.type_;
  }
  return *this;
}

// ###################################################################
//  Assignments
/**Assigns an arbitrary sequence of bytes value.*/
chi_data_types::Varying&
chi_data_types::Varying::operator=(const std::vector<std::byte>& value)
{
  type_ = VaryingDataType::ARBITRARY_BYTES;
  raw_data_ = value;
  data_initialized_ = true;
  return *this;
}

/**Assigns a string value.*/
chi_data_types::Varying&
chi_data_types::Varying::operator=(const std::string& value)
{
  type_ = VaryingDataType::STRING;
  PopulateRaw<std::string>(value);
  return *this;
}

// ###################################################################
//  Get values
/**Returns the string value if valid. Otherwise throws std::logic_error.*/
std::string chi_data_types::Varying::StringValue() const
{
  CheckTypeMatch(type_, VaryingDataType::STRING);
  CheckDataInitialized();

  if (reinterpret_cast<const char*>(raw_data_.data()) == nullptr)
    return std::string();
  return std::string(reinterpret_cast<const char*>(raw_data_.data()));
}

/**Returns the bool value if valid. Otherwise throws std::logic_error.*/
bool chi_data_types::Varying::BoolValue() const
{
  CheckTypeMatch(type_, VaryingDataType::BOOL);
  CheckDataInitialized();

  return *reinterpret_cast<const bool*>(&raw_data_[0]);
}

/**Returns the integer value if valid. Otherwise throws std::logic_error.*/
int64_t chi_data_types::Varying::IntegerValue() const
{
  CheckTypeMatch(type_, VaryingDataType::INTEGER);
  CheckDataInitialized();

  return *reinterpret_cast<const int64_t*>(&raw_data_[0]);
}

/**Returns the float value if valid. Otherwise throws std::logic_error.*/
double chi_data_types::Varying::FloatValue() const
{
  CheckTypeMatch(type_, VaryingDataType::FLOAT);
  CheckDataInitialized();

  return *reinterpret_cast<const double*>(&raw_data_[0]);
}

// ###################################################################
/**Returns the raw byte size associated with the type.*/
size_t chi_data_types::Varying::ByteSize() const { return raw_data_.size(); }

// ###################################################################
/**Returns a string value for the value.*/
std::string chi_data_types::Varying::PrintStr() const
{
  std::stringstream outstr;

  outstr << *this;

  return outstr.str();
}

// ###################################################################
/**Stream operator*/
std::ostream& operator<<(std::ostream& outstr,
                         const chi_data_types::Varying& value)
{
  if (value.Type() == chi_data_types::VaryingDataType::STRING)
    outstr << "\"" << value.StringValue() << "\"";
  else if (value.Type() == chi_data_types::VaryingDataType::FLOAT)
    outstr << value.FloatValue() << "(double)";
  else if (value.Type() == chi_data_types::VaryingDataType::INTEGER)
    outstr << value.IntegerValue();
  else if (value.Type() == chi_data_types::VaryingDataType::BOOL)
    outstr << (value.BoolValue() ? "true" : "false");
  return outstr;
}