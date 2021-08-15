#ifndef CHI_DATA_TYPES_VARYING_H
#define CHI_DATA_TYPES_VARYING_H

#include "chi_data_types.h"

#include <iostream>
#include <vector>
#include <memory>

namespace chi_data_types
{
class Varying
{
public:
  typedef char RawByte;
private:
  /**Byte size of the raw data*/
  size_t                m_num_bytes = 0;
  /**Raw byte-value data*/
  std::vector<RawByte>  m_raw_data;
  /**Flag indicating whether initialized or not*/
  bool                  m_data_initialized = false;
  /**Type specification*/
  VaryingDataType       m_type = VaryingDataType::VOID;

private:
  /**Utility that converts a type to a byte vector provided that
   * it has the sizeof() function defined for it.*/
  template<typename T> void PopulateRaw(const T& value)
    {
    auto src = reinterpret_cast<const RawByte*>(&value);

    m_num_bytes = sizeof(T);
    m_raw_data.resize(m_num_bytes);

    std::copy(src, src + m_num_bytes, &m_raw_data[0]);

    m_data_initialized = true;
    }

public:
  //Constructors
  /**Constructor for a string value.*/
  explicit Varying(const std::string& value);
  /**Constructor for a bool value.*/
  explicit Varying(const bool& value);
  /**Constructor for an integer value.*/
  explicit Varying(const int64_t& value);
  /**Constructor for a floating point value.*/
  explicit Varying(const double& value);

  /**Copy constructor.*/
  Varying(const Varying& other);

  /**Move constructor.*/
  Varying(Varying&& other) noexcept;

public:
  //Assignment operators
  /**Assignment operator. i.e., type_A = type_B*/
  Varying& operator=(const Varying& other);

  /**Assigns a string value.*/
  Varying& operator=(const std::string& value);
  /**Assigns a bool value.*/
  Varying& operator=(const bool& value);
  /**Assigns an integer value.*/
  Varying& operator=(const int64_t& value);
  /**Assign a floating point value.*/
  Varying& operator=(const double& value);

private:
  /**Checks if two VaryingDataType values match.
 * Type A is matched against type B.*/
  static void CheckTypeMatch(VaryingDataType type_A,
                             VaryingDataType type_B_required);
  /**Checks whether the data has been initialized.*/
  void CheckDataInitialized() const;

public:
  /**Returns the string value if valid. Otherwise throws std::logic_error.*/
  std::string StringValue() const;
  /**Returns the bool value if valid. Otherwise throws std::logic_error.*/
  bool        BoolValue() const;
  /**Returns the integer value if valid. Otherwise throws std::logic_error.*/
  int64_t     IntegerValue() const;
  /**Returns the float value if valid. Otherwise throws std::logic_error.*/
  double      FloatValue() const;

public:
  /**Returns the current-type of the variable.*/
  VaryingDataType Type() const {return m_type;}

public:
  ~Varying() = default;
};//class Varying
}//namespace chi_data_types

#endif //CHI_DATA_TYPES_VARYING_H