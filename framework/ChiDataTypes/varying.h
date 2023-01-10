#ifndef CHI_DATA_TYPES_VARYING_H
#define CHI_DATA_TYPES_VARYING_H

#include "chi_data_types.h"

#include <cstddef>
#include <iostream>
#include <vector>
#include <memory>

namespace chi_data_types
{
class Varying
{
private:
  /**Raw byte-value data*/
  std::vector<std::byte>  m_raw_data;
  /**Flag indicating whether initialized or not*/
  bool                    m_data_initialized = false;
  /**Type specification*/
  VaryingDataType         m_type      = VaryingDataType::VOID;
  std::string             m_type_name = "VOID";

private:
  /**Utility that converts a type to a byte vector provided that
   * it has the sizeof() function defined for it.*/
  template<typename T> void PopulateRaw(const T& value)
  {
    auto src = reinterpret_cast<const std::byte*>(&value);

    size_t num_bytes = sizeof(T);
    m_raw_data.resize(num_bytes);

    std::copy(src, src + num_bytes, &m_raw_data[0]);

    m_data_initialized = true;
  }
  /**Checks if two VaryingDataType values match.
   * Type A is matched against type B.*/
  static void CheckTypeMatch(VaryingDataType type_A,
                             VaryingDataType type_B_required);
  /**Checks whether the data has been initialized.*/
  void CheckDataInitialized() const;

public:
  //Constructors
  /**Constructor for an arbitrary sequence of bytes value.*/
  explicit Varying(const std::vector<std::byte>& value);
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
  //Copy assignment operator
  /**Assignment operator. i.e., type_A = type_B*/
  Varying& operator=(const Varying& other);

  //Assignment operators
  /**Assigns an arbitrary sequence of bytes value.*/
  Varying& operator=(const std::vector<std::byte>& value);
  /**Assigns a string value.*/
  Varying& operator=(const std::string& value);
  /**Assigns a bool value.*/
  Varying& operator=(const bool& value);
  /**Assigns an integer value.*/
  Varying& operator=(const int64_t& value);
  /**Assign a floating point value.*/
  Varying& operator=(const double& value);


public:
  /**Returns the string value if valid. Otherwise throws std::logic_error.*/
  std::string StringValue() const;
  /**Returns the bool value if valid. Otherwise throws std::logic_error.*/
  bool        BoolValue() const;
  /**Returns the integer value if valid. Otherwise throws std::logic_error.*/
  int64_t     IntegerValue() const;
  /**Returns the float value if valid. Otherwise throws std::logic_error.*/
  double      FloatValue() const;

  /**Returns the raw byte size associated with the type.*/
  size_t      ByteSize() const;


public:
  /**Returns the current-type of the variable.*/
  VaryingDataType Type() const {return m_type;}
  /**Returns the string type name of the type.*/
  std::string TypeName() const {return m_type_name;}

public:
  ~Varying() = default;
};//class Varying
}//namespace chi_data_types

#endif //CHI_DATA_TYPES_VARYING_H