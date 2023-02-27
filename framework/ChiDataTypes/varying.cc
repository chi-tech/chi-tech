#include "varying.h"

#include <algorithm>

/**Provides a string-name for an enumerated VaryingDataType.*/
std::string chi_data_types::
  VaryingDataTypeStringName(chi_data_types::VaryingDataType type)
{
  switch (type)
  {
    case VaryingDataType::VOID:            return "VOID";
    case VaryingDataType::ARBITRARY_BYTES: return "ARBITRARY_BYTES";
    case VaryingDataType::STRING:          return "STRING";
    case VaryingDataType::BOOL:            return "BOOL";
    case VaryingDataType::INTEGER:         return "INTEGER";
    case VaryingDataType::FLOAT:           return "FLOAT";
    default: return "UNKNOWN";
  }
}

//###################################################################
/**PopulateRaw template specialization for std::string.*/
template<>
void chi_data_types::Varying::PopulateRaw(const std::string& value)
{
  raw_data_.resize(value.size());
  std::transform(value.begin(), value.end(), raw_data_.begin(),
                 [](char c){return std::byte(c);});
  if (value.back() != '\0')
    raw_data_.push_back(std::byte('\0'));

  data_initialized_ = true;
}

//################################################################### Constructors
/**Constructor for an arbitrary sequence of bytes value.*/
chi_data_types::Varying::
  Varying(const std::vector<std::byte>& value) :
  type_(VaryingDataType::ARBITRARY_BYTES),
  type_name_(VaryingDataTypeStringName(type_))
{
  raw_data_ = value;
  data_initialized_ = true;
}

/**Constructor for a string value.*/
chi_data_types::Varying::
  Varying(const std::string& value) :
  type_(VaryingDataType::STRING),
  type_name_(VaryingDataTypeStringName(type_))
{ PopulateRaw<std::string>(value); }

/**Constructor for a bool value.*/
chi_data_types::Varying::
  Varying(const bool& value) :
  type_(VaryingDataType::BOOL),
  type_name_(VaryingDataTypeStringName(type_))
{ PopulateRaw<bool>(value); }

/**Constructor for an integer value.*/
chi_data_types::Varying::
  Varying(const int64_t& value) :
  type_(VaryingDataType::INTEGER),
  type_name_(VaryingDataTypeStringName(type_))
{ PopulateRaw<int64_t>(value); }

/**Constructor for a floating point value.*/
chi_data_types::Varying::
  Varying(const double& value) :
  type_(VaryingDataType::FLOAT),
  type_name_(VaryingDataTypeStringName(type_))
{ PopulateRaw<double>(value); }

/**Copy constructor.*/
chi_data_types::Varying::Varying(const Varying& other)
{
  raw_data_         = other.raw_data_;
  data_initialized_ = other.data_initialized_;
  type_             = other.type_;
  type_name_        = other.type_name_;
}

/**Move constructor.*/
chi_data_types::Varying::Varying(Varying&& other) noexcept
{
  std::swap(raw_data_, other.raw_data_);
  std::swap(data_initialized_, other.data_initialized_);
  std::swap(type_, other.type_);
  std::swap(type_name_, other.type_name_);
}

/**Assignment operator. i.e., type_A = type_B*/
chi_data_types::Varying& chi_data_types::Varying::operator=(const Varying& other)
{
  if (this != &other)
  {
    raw_data_         = other.raw_data_;
    data_initialized_ = other.data_initialized_;
    type_             = other.type_;
    type_name_        = other.type_name_;
  }
  return *this;
}

//################################################################### Assignments
/**Assigns an arbitrary sequence of bytes value.*/
chi_data_types::Varying&
  chi_data_types::Varying::operator=(const std::vector<std::byte> &value)
{
  type_ = VaryingDataType::ARBITRARY_BYTES;
  type_name_ = VaryingDataTypeStringName(type_);
  raw_data_ = value;
  data_initialized_ = true;
  return *this;
}

/**Assigns a string value.*/
chi_data_types::Varying&
  chi_data_types::Varying::operator=(const std::string& value)
{
  type_ = VaryingDataType::STRING;
  type_name_ = VaryingDataTypeStringName(type_);
  PopulateRaw<std::string>(value);
  return *this;
}

/**Assigns a bool value.*/
chi_data_types::Varying& chi_data_types::Varying::operator=(const bool& value)
{
  type_ = VaryingDataType::BOOL;
  type_name_ = VaryingDataTypeStringName(type_);
  PopulateRaw<bool>(value);

  return *this;
}

/**Assigns an integer value.*/
chi_data_types::Varying& chi_data_types::Varying::operator=(const int64_t& value)
{
  type_ = VaryingDataType::INTEGER;
  type_name_ = VaryingDataTypeStringName(type_);
  PopulateRaw<int64_t>(value);

  return *this;
}

/**Assign a floating point value.*/
chi_data_types::Varying& chi_data_types::Varying::operator=(const double& value)
{
  type_ = VaryingDataType::FLOAT;
  type_name_ = VaryingDataTypeStringName(type_);
  PopulateRaw<double>(value);

  return *this;
}

//################################################################### Get values

/**Checks if two VaryingDataType values match.
 * Type A is matched against type B.*/
void chi_data_types::Varying::
  CheckTypeMatch(const VaryingDataType type_A,
                 const VaryingDataType type_B_required)
{
  if (type_A != type_B_required)
    throw std::logic_error("Varying data type does not correspond to "
                           "the required type," +
                           VaryingDataTypeStringName(type_B_required));
}

/**Checks whether the data has been initialized.*/
void chi_data_types::Varying::CheckDataInitialized() const
{
  if (not data_initialized_)
    throw std::logic_error("Varying data type used uninitialized.");
}

/**Returns the string value if valid. Otherwise throws std::logic_error.*/
std::string chi_data_types::Varying::StringValue() const
{
  CheckTypeMatch(type_, VaryingDataType::STRING);
  CheckDataInitialized();

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

/**Returns the raw byte size associated with the type.*/
size_t chi_data_types::Varying::ByteSize() const
{
  return raw_data_.size();
}