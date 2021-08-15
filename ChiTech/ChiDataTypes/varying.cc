#include "varying.h"

/**Provides a string-name for an enumerated VaryingDataType.*/
std::string chi_data_types::
  VaryingDataTypeStringName(chi_data_types::VaryingDataType type)
{
  switch (type)
  {
    case VaryingDataType::VOID:    return "VOID";
    case VaryingDataType::STRING:  return "STRING";
    case VaryingDataType::BOOL:    return "BOOL";
    case VaryingDataType::INTEGER: return "INTEGER";
    case VaryingDataType::FLOAT:   return "FLOAT";
    default: return "UNKNOWN";
  }
}

//################################################################### Constructors
/**Constructor for a string value.*/
chi_data_types::Varying::Varying(const std::string& value) :
  m_data_initialized(true),
  m_type(VaryingDataType::STRING)
{
  m_raw_data  = std::vector<RawByte>(value.begin(), value.end());
  m_raw_data.push_back('\0'); //Necessary to complete c_str()
  m_num_bytes = m_raw_data.size();
}

/**Constructor for a bool value.*/
chi_data_types::Varying::Varying(const bool& value) :
  m_type(VaryingDataType::BOOL)
{
  PopulateRaw<bool>(value);
}

/**Constructor for an integer value.*/
chi_data_types::Varying::Varying(const int64_t& value) :
  m_type(VaryingDataType::INTEGER)
{
  PopulateRaw<int64_t>(value);
}

/**Constructor for a floating point value.*/
chi_data_types::Varying::Varying(const double& value) :
  m_type(VaryingDataType::FLOAT)
{
  PopulateRaw<double>(value);
}

/**Copy constructor.*/
chi_data_types::Varying::Varying(const Varying& other)
{
  m_type             = other.m_type;
  m_num_bytes        = other.m_num_bytes;
  m_raw_data         = other.m_raw_data;
  m_data_initialized = other.m_data_initialized;
}

/**Move constructor.*/
chi_data_types::Varying::Varying(Varying&& other) noexcept
{
  std::swap(m_type, other.m_type);
  std::swap(m_num_bytes, other.m_num_bytes);
  std::swap(m_raw_data, other.m_raw_data);
  std::swap(m_data_initialized, other.m_data_initialized);
}

/**Assignment operator. i.e., type_A = type_B*/
chi_data_types::Varying& chi_data_types::Varying::operator=(const Varying& other)
{
  if (this != &other)
  {
    m_type             = other.m_type;
    m_num_bytes        = other.m_num_bytes;
    m_raw_data         = other.m_raw_data;
    m_data_initialized = other.m_data_initialized;
  }
  return *this;
}

//################################################################### Assignments

/**Assigns a string value.*/
chi_data_types::Varying&
  chi_data_types::Varying::operator=(const std::string& value)
{
  m_type = VaryingDataType::STRING;

  m_raw_data  = std::vector<RawByte>(value.begin(), value.end());
  m_raw_data.push_back('\0');
  m_num_bytes = m_raw_data.size();

  m_data_initialized = true;
  return *this;
}

/**Assigns a bool value.*/
chi_data_types::Varying& chi_data_types::Varying::operator=(const bool& value)
{
  m_type = VaryingDataType::BOOL;
  PopulateRaw<bool>(value);

  return *this;
}

/**Assigns an integer value.*/
chi_data_types::Varying& chi_data_types::Varying::operator=(const int64_t& value)
{
  m_type = VaryingDataType::INTEGER;
  PopulateRaw<int64_t>(value);

  return *this;
}

/**Assign a floating point value.*/
chi_data_types::Varying& chi_data_types::Varying::operator=(const double& value)
{
  m_type = VaryingDataType::FLOAT;
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
  if (not m_data_initialized)
    throw std::logic_error("Varying data type used uninitialized.");
}

/**Returns the string value if valid. Otherwise throws std::logic_error.*/
std::string chi_data_types::Varying::StringValue() const
{
  CheckTypeMatch(m_type, VaryingDataType::STRING);
  CheckDataInitialized();

  return std::string(m_raw_data.data());
}

/**Returns the bool value if valid. Otherwise throws std::logic_error.*/
bool chi_data_types::Varying::BoolValue() const
{
  CheckTypeMatch(m_type, VaryingDataType::BOOL);
  CheckDataInitialized();

  return *reinterpret_cast<const bool*>(&m_raw_data[0]);
}

/**Returns the integer value if valid. Otherwise throws std::logic_error.*/
int64_t chi_data_types::Varying::IntegerValue() const
{
  CheckTypeMatch(m_type, VaryingDataType::INTEGER);
  CheckDataInitialized();

  return *reinterpret_cast<const int64_t*>(&m_raw_data[0]);
}

/**Returns the float value if valid. Otherwise throws std::logic_error.*/
double chi_data_types::Varying::FloatValue() const
{
  CheckTypeMatch(m_type, VaryingDataType::FLOAT);
  CheckDataInitialized();

  return *reinterpret_cast<const double*>(&m_raw_data[0]);
}