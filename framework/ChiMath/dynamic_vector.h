#ifndef _chi_math_DVectorNX_h
#define _chi_math_DVectorNX_h

#include <vector>
#include <stdexcept>
#include <sstream>

namespace chi_mesh
{
  struct Vector3;
}

namespace chi_math
{
  template<class NumberFormat>
  class DynamicVector;
}

//###################################################################
/** General dynamic vector utility.*/
template<class NumberFormat>
class chi_math::DynamicVector
{
public:
  std::vector<NumberFormat> elements;

  /**Default constructor. Does nothing.*/
  DynamicVector()
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are "
                  "supported for DVectorNX." );
  }

  /**Constructor with number of entries. Value defaults.*/
  DynamicVector(size_t N) :
    elements(N)
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are "
                  "supported for DVectorNX." );
  }

  /**Constructor with number of entries and default value.*/
  DynamicVector(size_t N, NumberFormat value) :
    elements(N, value)
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are "
                  "supported for DVectorNX." );
  }

  /**Copy constructor.*/
  DynamicVector(const DynamicVector& other) { elements = other.elements;}

  /**Assignment operator.*/
  DynamicVector& operator=(const DynamicVector& other)
  {
    elements = other.elements;
    return *this;
  }

  /**Move constructor.*/
  DynamicVector(DynamicVector&& other) { elements = std::move(other.elements);}

  /**Move assignment operator.*/
  DynamicVector& operator=(DynamicVector&& other)
  {
    elements = std::move(other.elements);
    return *this;
  }

  /**Constructor with vector.*/
  DynamicVector(const std::vector<double>& in) { elements = in;}

  /**Copy constructor with vector.*/
  DynamicVector& operator=(const std::vector<double>& in)
  {
    elements = in;
    return *this;
  }

  /**Constructor with vector.*/
  DynamicVector(std::initializer_list<NumberFormat> in) { elements = in;}

  /**Copy constructor with vector.*/
  DynamicVector& operator=(std::initializer_list<NumberFormat> in)
  {
    elements = in;
    return *this;
  }

  //============================================= Element access
  NumberFormat& operator[](size_t i) {return elements[i];}
  const NumberFormat& operator[](size_t i) const {return elements[i];}

  NumberFormat& at(size_t i) {return elements.at(i);}

  NumberFormat& back() {return elements.back();}

  NumberFormat& front() {return elements.front();}

  NumberFormat* data() {return elements.data();}

  void clear() {elements.clear();}

  void resize(size_t dim) {elements.resize(dim); }
  void resize(size_t dim, const NumberFormat& val)
  {
    elements.resize(dim, val);
    for (auto& entry : elements)
      entry = val;
  }

  void reserve(size_t dim) {elements.reserve(dim);}

  void push_back(const NumberFormat& val) {elements.push_back(val);}
  void pop_back() {elements.pop_back();}

  bool empty() const noexcept {return elements.empty();}

  //============================================= Iterator access
  typename std::vector<NumberFormat>::iterator begin() {return elements.begin();}

  typename std::vector<NumberFormat>::iterator end() {return elements.end();}

  size_t size() const {return elements.size();}

  void bounds_check(const size_t a, const size_t b) const
  {
    if (a != b)
    {
      std::length_error excp("Mismatched sizes of DVectorNX");
      throw excp;
    }
  }

  //============================================= Addition
  /**Component-wise addition of two vectors.
   * \f$ \vec{w} = \vec{x} + \vec{y} \f$*/
  DynamicVector operator+(const DynamicVector& rhs) const
  {
    bounds_check(size(),rhs.size());
    DynamicVector<NumberFormat> newVector(size());
    for (int i = 0; i<size();++i)
      newVector.elements[i] = elements[i] + rhs.elements[i];

    return newVector;
  }

  /**In-place component-wise addition of two vectors.
 * \f$ \vec{x} = \vec{x} + \vec{y} \f$*/
  DynamicVector& operator+=(const DynamicVector& rhs)
  {
    bounds_check(size(),rhs.size());
    for (int i = 0; i<size();++i)
      elements[i] += rhs.elements[i];

    return *this;
  }

  //=========================================== Subtraction
  /**Component-wise subtraction.
 * \f$ \vec{w} = \vec{x} - \vec{y} \f$*/
  DynamicVector operator-(const DynamicVector& rhs) const
  {
    bounds_check(size(),rhs.size());
    DynamicVector<NumberFormat> newVector(size());
    for (int i = 0; i<size();++i)
      newVector.elements[i] = elements[i] - rhs.elements[i];

    return newVector;
  }

  /**In-place component-wise subtraction.
 * \f$ \vec{x} = \vec{x} - \vec{y} \f$*/
  DynamicVector& operator-=(const DynamicVector& rhs)
  {
    bounds_check(size(),rhs.size());
    for (int i = 0; i<size();++i)
      elements[i] -= rhs.elements[i];

    return *this;
  }

  //=========================================== Multiplication
  /**Vector component-wise multiplication by scalar.
 * \f$ \vec{w} = \vec{x} \alpha \f$*/
  DynamicVector operator*(const NumberFormat value) const
  {
    DynamicVector<NumberFormat> newVector(size());
    for (int i = 0;i<size();++i)
      newVector.elements[i] = elements[i] * value;

    return newVector;
  }

  /**Vector in-place component-wise multiplication by scalar.
 * \f$ \vec{x} = \vec{x} \alpha \f$*/
  DynamicVector& operator*=(const NumberFormat value)
  {
    for (int i = 0;i<size();++i)
      elements[i] *= value;

    return *this;
  }

  /**Vector component-wise multiplication.
 * \f$ w_i = x_i y_i \f$*/
  DynamicVector operator*(const DynamicVector& rhs) const
  {
    bounds_check(size(),rhs.size());
    DynamicVector<NumberFormat> newVector(size());
    for (int i = 0; i<size(); ++i)
      newVector.elements[i] = elements[i] * rhs.elements[i];

    return newVector;
  }

  /**Vector in-place component-wise multiplication.
 * \f$ x_i = x_i y_i \f$*/
  DynamicVector& operator*=(const DynamicVector& rhs)
  {
    bounds_check(size(),rhs.size());
    for (int i = 0; i<size(); ++i)
      elements[i] *= rhs.elements[i];

    return *this;
  }

  //=========================================== Division
  /**Vector component-wise division by scalar.
 * \f$ w_i = \frac{x_i}{\alpha} \f$*/
  DynamicVector operator/(const NumberFormat value) const
  {
    DynamicVector<NumberFormat> newVector(size());
    for (int i = 0;i<size();++i)
      newVector.elements[i] = elements[i] / value;

    return newVector;
  }

  /**Vector in-place component-wise division by scalar.
 * \f$ x_i = \frac{x_i}{\alpha} \f$*/
  DynamicVector& operator/=(const NumberFormat value)
  {
    for (int i = 0;i<size();++i)
      elements[i] /= value;

    return *this;
  }

  /**Vector component-wise division.
 * \f$ w_i = \frac{x_i}{y_i} \f$*/
  DynamicVector operator/(const DynamicVector& rhs) const
  {
    bounds_check(size(),rhs.size());
    DynamicVector<NumberFormat> newVector(size());
    for (int i = 0; i<size(); ++i)
      newVector.elements[i] = elements[i] / rhs.elements[i];

    return newVector;
  }

  /**Vector in-place component-wise division.
  * \f$ x_i = \frac{x_i}{y_i} \f$*/
  DynamicVector& operator/=(const DynamicVector& rhs)
  {
    bounds_check(size(),rhs.size());
    for (int i = 0; i<size(); ++i)
      elements[i] /= rhs.elements[i];

    return *this;
  }

  //============================================= Operations
  /**Vector dot-product.
   * \f$ \vec{w} = \vec{x} \bullet \vec{y} \f$ */
  NumberFormat Dot(const DynamicVector& rhs) const
  {
    bounds_check(size(),rhs.size());
    NumberFormat value = 0.0;
    for (int i = 0; i<size(); ++i)
      value += elements[i]*rhs.elements[i];

    return value;
  }

  /**Computes the L2-norm of the vector. Otherwise known as the length of
    * a 3D vector.*/
  NumberFormat Norm() const
  {
    NumberFormat value = 0.0;
    for (int i = 0; i<size();++i)
      value += elements[i]*elements[i];

    value = sqrt(value);
    return value;
  }

  /**Computes the square of the L2-norm of the vector. This eliminates the
  * usage of the square root and is therefore less expensive that a proper
  * L2-norm. Useful if only comparing distances.*/
  NumberFormat NormSquare() const
  {
    NumberFormat value = 0.0;
    for (int i = 0; i<size();++i)
      value += elements[i]*elements[i];

    return value;
  }

  /**Normalizes the vector in-place.*/
  void Normalize()
  {
    NumberFormat norm = this->Norm();
    for (int i = 0;i<size();++i)
      elements[i] = elements[i]/norm;
  }

  /**Returns a normalized version of the vector.*/
  DynamicVector Normalized() const
  {
    NumberFormat norm = this->Norm();
    DynamicVector<NumberFormat> newVector(size());
    for (int i = 0;i<size();++i)
      newVector.elements[i] = elements[i]/norm;

    return newVector;
  }

  /**Prints the vector to a string and then returns the string.*/
  std::string PrintStr() const
  {
    std::stringstream out;
    out<<"[";
    for (int i = 0; i<(size()-1); ++i)
      out<<elements[i]<<" ";
    out<<elements[size()-1]<<"]";

    return out.str();
  }
};

/**Multiplication by a scalar from the left.*/
template<class NumberFormat>
chi_math::DynamicVector<NumberFormat>
operator*(const double value,const chi_math::DynamicVector<NumberFormat>& that)
{
  chi_math::DynamicVector<NumberFormat> newVector(that.size());
  for (int i = 0; i<that.size();++i)
    newVector.elements[i] = that.elements[i]*value;
  return newVector;
}

#endif