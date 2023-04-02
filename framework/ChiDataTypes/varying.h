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
  std::vector<std::byte>  raw_data_;
  /**Flag indicating whether initialized or not*/
  bool                    data_initialized_ = false;
  /**Type specification*/
  VaryingDataType         type_      = VaryingDataType::VOID;

private:
  /**Utility that converts a type to a byte vector provided that
   * it has the sizeof() function defined for it.*/
  template<typename T>
  void PopulateRaw(const T& value)
  {
    auto src = reinterpret_cast<const std::byte*>(&value);

    size_t num_bytes = sizeof(T);
    raw_data_.resize(num_bytes);

    std::copy(src, src + num_bytes, &raw_data_[0]);

    data_initialized_ = true;
  }

  /**Checks if two VaryingDataType values match.
   * Type A is matched against type B.*/
  void CheckTypeMatch(VaryingDataType type_A,
                      VaryingDataType type_B_required) const;
  /**Checks whether the data has been initialized.*/
  void CheckDataInitialized() const;

public:
  //Helpers
  template<typename T> struct IsBool
  {static constexpr bool value = std::is_same_v<T,bool>;};
  template<typename T> struct IsFloat
  {static constexpr bool value = std::is_floating_point_v<T>;};
  template<typename T> struct IsInteger
  {static constexpr bool value = std::is_integral_v<T> and
                                 not std::is_same_v<T,bool>;};

  template<typename T> using BoolType =
    typename std::enable_if_t<IsBool<T>::value, T>;
  template<typename T> using FloatType =
    typename std::enable_if_t<IsFloat<T>::value, T>;
  template<typename T> using IntegerType =
    typename std::enable_if_t<IsInteger<T>::value, T>;

private:
  template<typename T>
  BoolType<T> CastValue(const T& value) { return value; }

  template<typename T>
  FloatType<T> CastValue(const T& value) { return static_cast<double>(value);}

  template<typename T>
  IntegerType<T> CastValue(const T& value) {return static_cast<int64_t>(value);}

public:
  /**Generalized constructor for bool, integral- and float-types. This
   * constructor has been specialized for std::string and
   * std::vector<std::byte>.*/
  template<typename T>
    explicit Varying(const T& value)
  {
    constexpr bool is_supported_type =
      IsBool<T>::value or IsFloat<T>::value or IsInteger<T>::value;
    static_assert(is_supported_type,
      "Constructor called with unsupported type");

    if (IsBool<T>::value) type_ = VaryingDataType::BOOL;
    else if (IsFloat<T>::value) type_ = VaryingDataType::FLOAT;
    else if (IsInteger<T>::value) type_ = VaryingDataType::INTEGER;

    PopulateRaw<T>(CastValue(value));
  }

  //Constructors
  /**Constructor for an arbitrary sequence of bytes value.*/
  explicit Varying(const std::vector<std::byte>& value);
  /**Constructor for a string value.*/
  explicit Varying(const std::string& value);

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
  template<typename T, std::enable_if_t<IsBool<T>::value, bool> = true>
  Varying& operator=(const T& value)
  {
    type_ = VaryingDataType::BOOL;
    PopulateRaw<bool>(value);

    return *this;
  }

  /**Assigns an integer value.*/
  template<typename T, std::enable_if_t<IsInteger<T>::value, bool> = true>
  Varying& operator=(const T& value)
  {
    type_ = VaryingDataType::INTEGER;
    PopulateRaw<int64_t>(static_cast<int64_t>(value));

    return *this;
  }

  /**Assign a floating point value.*/
  template<typename T, std::enable_if_t<IsFloat<T>::value, bool> = true>
  Varying& operator=(const T& value)
  {
    type_ = VaryingDataType::FLOAT;
    PopulateRaw<double>(static_cast<double>(value));

    return *this;
  }


public:
  //More Helpers
  template<typename T> struct IsString
  {static constexpr bool value = std::is_same_v<T,std::string>;};
  template<typename T> struct IsSignedInteger
  {static constexpr bool value = std::is_integral_v<T> and
                                 std::is_signed_v<T> and
                                 not std::is_same_v<T,bool>;};
  template<typename T> struct IsUnsignedInteger
  {static constexpr bool value = std::is_integral_v<T> and
                                 std::is_unsigned_v<T> and
                                 not std::is_same_v<T,bool>;};

  template<typename T> using StringType =
    typename std::enable_if_t<IsString<T>::value, T>;
  template<typename T> using SignedIntegerType =
    typename std::enable_if_t<IsSignedInteger<T>::value, T>;
  template<typename T> using UnsignedIntegerType =
    typename std::enable_if_t<IsUnsignedInteger<T>::value, T>;

  /**Returns values of type bool if able.*/
  template<typename T>
  BoolType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::BOOL);
    CheckDataInitialized();

    return *reinterpret_cast<const bool*>(&raw_data_[0]);
  }

  /**Returns floating point values if able.*/
  template<typename T>
  FloatType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::FLOAT);
    CheckDataInitialized();

    const double value = *reinterpret_cast<const double*>(&raw_data_[0]);

    return static_cast<T>(value);
  }

  /**Returns a string if able.*/
  template<typename T>
  StringType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::STRING);
    CheckDataInitialized();

    return std::string(reinterpret_cast<const char*>(raw_data_.data()));
  }

  /**Returns a signed integer if able.*/
  template<typename T>
  SignedIntegerType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::INTEGER);
    CheckDataInitialized();

    const int64_t value = *reinterpret_cast<const int64_t*>(&raw_data_[0]);

    return static_cast<T>(value);
  }

  /**Returns an unsigned integer if able.*/
  template<typename T>
  UnsignedIntegerType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::INTEGER);
    CheckDataInitialized();

    const int64_t value = *reinterpret_cast<const int64_t*>(&raw_data_[0]);

    if (value < 0)
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                             ": Attempt to cast negative number to unsigned.");

    return static_cast<T>(value);
  }

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
  VaryingDataType Type() const {return type_;}
  /**Returns the string type name of the type.*/
  std::string TypeName() const {return VaryingDataTypeStringName(type_);}

public:
  ~Varying() = default;
};//class Varying
}//namespace chi_data_types

#endif //CHI_DATA_TYPES_VARYING_H