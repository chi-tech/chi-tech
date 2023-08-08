#ifndef CHITECH_ALLOWABLE_RANGE_H
#define CHITECH_ALLOWABLE_RANGE_H

#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <sstream>

namespace chi_data_types
{

// ##################################################################
/**Base class for an allowable range.*/
class AllowableRange
{
protected:
  virtual bool ChildIsAllowable(Varying value) const = 0;
  virtual std::string AllowableRangeStr() const = 0;

public:
  /**Returns `true` is the value is within the allowable range.*/
  template <typename T>
  bool IsAllowable(const T& value)
  {
    return ChildIsAllowable(Varying(value));
  }

  /**Provides a formatted string of the message to display if a
   * value is out of range.*/
  template <typename T>
  std::string OutOfRangeString(const std::string& scope, const T& value) const
  {
    const Varying vvalue(value);
    std::stringstream outstr;

    if (not scope.empty()) outstr << "Parameter \"" << scope << "\": ";

    outstr << "Value " << vvalue << " out of range. ";
    outstr << "Constraints: " << AllowableRangeStr() << ".";

    return outstr.str();
  }

  /**Prints the allowable constraints to a string.*/
  std::string PrintRange() { return AllowableRangeStr(); }

  virtual ~AllowableRange() = default;
};

// ##################################################################
/**Range comprising a list of values.*/
class AllowableRangeList : public AllowableRange
{
protected:
  std::vector<Varying> list_;

public:
  template <typename T>
  AllowableRangeList(const std::initializer_list<T>& raw_list)
  {
    for (const auto& val : raw_list)
      list_.emplace_back(val);
  }

  template <typename T>
  AllowableRangeList(const std::vector<T>& raw_list)
  {
    for (const auto& val : raw_list)
      list_.emplace_back(val);
  }

  template <typename T>
  static std::unique_ptr<AllowableRangeList>
  New(const std::initializer_list<T>& raw_list)
  {
    return std::unique_ptr<AllowableRangeList>{
      new AllowableRangeList(raw_list)};
  }

  template <typename T>
  static std::unique_ptr<AllowableRangeList>
  New(const std::vector<T>& raw_list)
  {
    return std::unique_ptr<AllowableRangeList>{
      new AllowableRangeList(raw_list)};
  }

protected:
  bool ChildIsAllowable(Varying value) const override
  {
    return std::find(list_.begin(), list_.end(), value) != list_.end();
  }
  std::string AllowableRangeStr() const override
  {
    std::stringstream outstr;
    for (const auto& value : list_)
    {
      outstr << value;
      if (value != list_.back()) outstr << ", ";
    }
    return outstr.str();
  }
};

class AllowableRangeLowHighLimit;

// ##################################################################
/**Lower limit range.*/
class AllowableRangeLowLimit : public AllowableRange
{
protected:
  const Varying low_limit_;
  bool low_closed_ = true;

public:
  template <typename T>
  explicit AllowableRangeLowLimit(const T& low_value, bool low_closed = true)
    : low_limit_(Varying(low_value)), low_closed_(low_closed)
  {
  }

  template <typename T>
  static std::unique_ptr<AllowableRangeLowLimit> New(const T& low_value,
                                                     bool low_closed = true)
  {
    return std::unique_ptr<AllowableRangeLowLimit>{
      new AllowableRangeLowLimit(low_value, low_closed)};
  }

protected:
  friend class AllowableRangeLowHighLimit;
  bool ChildIsAllowable(Varying value) const override
  {
    if (value.Type() != low_limit_.Type()) return false;

    if (low_closed_) return value >= low_limit_;
    else
      return value > low_limit_;
  }
  std::string AllowableRangeStr() const override
  {
    if (low_closed_) return std::string(">= ") + low_limit_.PrintStr();
    else
      return std::string("> ") + low_limit_.PrintStr();
  }
};

// ##################################################################
/**Upper limit range.*/
class AllowableRangeHighLimit : public AllowableRange
{
protected:
  const Varying hi_limit_;
  bool hi_closed_ = true;

public:
  template <typename T>
  explicit AllowableRangeHighLimit(const T& hi_value, bool hi_closed = true)
    : hi_limit_(Varying(hi_value)), hi_closed_(hi_closed)
  {
  }

  template <typename T>
  static std::unique_ptr<AllowableRangeHighLimit> New(const T& hi_value,
                                                      bool hi_closed = true)
  {
    return std::unique_ptr<AllowableRangeHighLimit>{
      new AllowableRangeHighLimit(hi_value, hi_closed)};
  }

protected:
  friend class AllowableRangeLowHighLimit;
  bool ChildIsAllowable(Varying value) const override
  {
    if (value.Type() != hi_limit_.Type()) return false;

    if (hi_closed_) return value <= hi_limit_;
    else
      return value < hi_limit_;
  }
  std::string AllowableRangeStr() const override
  {
    if (hi_closed_) return std::string("<= ") + hi_limit_.PrintStr();
    else
      return std::string("< ") + hi_limit_.PrintStr();
  }
};

// ##################################################################
/**Upper and lower limit range.*/
class AllowableRangeLowHighLimit : public AllowableRange
{
protected:
  AllowableRangeLowLimit low_range_;
  AllowableRangeHighLimit hi_range;

public:
  template <typename T>
  explicit AllowableRangeLowHighLimit(const T& low_value,
                                      const T& hi_value,
                                      bool low_closed = true,
                                      bool hi_closed = true)
    : low_range_(Varying(low_value), low_closed),
      hi_range(Varying(hi_value), hi_closed)
  {
  }

  template <typename T>
  static std::unique_ptr<AllowableRangeLowHighLimit> New(const T& low_value,
                                                         const T& hi_value,
                                                         bool low_closed = true,
                                                         bool hi_closed = true)
  {
    return std::unique_ptr<AllowableRangeLowHighLimit>{
      new AllowableRangeLowHighLimit(
        low_value, hi_value, low_closed, hi_closed)};
  }

protected:
  bool ChildIsAllowable(Varying value) const override
  {
    return low_range_.ChildIsAllowable(value) and
           hi_range.ChildIsAllowable(value);
  }
  std::string AllowableRangeStr() const override
  {

    return low_range_.AllowableRangeStr() + ", " + hi_range.AllowableRangeStr();
  }
};

} // namespace chi_data_types

#endif // CHITECH_ALLOWABLE_RANGE_H
