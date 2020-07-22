#ifndef _chi_math_dynamic_matrixNX_h
#define _chi_math_dynamic_matrixNX_h

#include <vector>
#include <stdexcept>

namespace chi_math
{
  template<class NumberFormat>
  class DMatrixNX;
}

#include "chi_math_dynamic_vectorNX.h"

//###################################################################
/**General dynamic matrix utility.*/
template<class NumberFormat>
class chi_math::DMatrixNX
{
  typedef std::pair<size_t,size_t> MatDim;
public:
  std::vector<std::vector<NumberFormat>> elements;

  /**Default constructor. Does nothing.*/
  DMatrixNX()
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are "
                  "supported for DMatrixNX." );
  }

  /**Constructor with number of entries. Value defaults.*/
  DMatrixNX(size_t Nrows,size_t Ncols) :
    elements(Nrows,std::vector<NumberFormat>(Ncols))
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are "
                  "supported for DMatrixNX." );
  }

  /**Constructor with number of entries and default value.*/
  DMatrixNX(size_t Nrows, size_t Ncols, NumberFormat value) :
    elements(Nrows, std::vector<NumberFormat>(Ncols,value))
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are "
                  "supported for DMatrixNX." );
  }

  /**Copy constructor.*/
  DMatrixNX(const DMatrixNX& other) { elements = other.elements;}

  /**Assignment operator.*/
  DMatrixNX& operator=(const DMatrixNX& other)
  {
    elements = other.elements;
    return *this;
  }

  /**Move constructor.*/
  DMatrixNX(DMatrixNX&& other) { elements = std::move(other.elements);}

  /**Move assignment operator.*/
  DMatrixNX& operator=(DMatrixNX&& other)
  {
    elements = std::move(other.elements);
    return *this;
  }

  /**Constructor with vector.*/
  DMatrixNX(const std::vector<std::vector<double>>& in) { elements = in;}

  /**Copy constructor with vector.*/
  DMatrixNX& operator=(const std::vector<std::vector<double>>& in)
  {
    elements = in;
    return *this;
  }

  /**Constructor with vector.*/
  DMatrixNX(std::initializer_list<std::initializer_list<NumberFormat>> in)
  {
    elements.clear();
    for (auto& v : in)
      elements.push_back(v);
  }

  /**Copy constructor with vector.*/
  DMatrixNX& operator=(std::initializer_list<std::initializer_list<NumberFormat>> in)
  {
    elements.clear();
    for (auto& v : in)
      elements.push_back(v);
    return *this;
  }

  //============================================= Element access
  std::vector<NumberFormat>& operator[](size_t i) {return elements[i];}

  std::vector<NumberFormat>& at(size_t i) {return elements.at(i);}

  std::vector<NumberFormat>& back() {return elements.back();}

  std::vector<NumberFormat>& front() {return elements.front();}

  std::vector<NumberFormat>* data() {return elements.data();}

  void clear() {elements.clear();}

  void resize(size_t Nrows, size_t Ncols)
  {
    elements.resize(Nrows);
    for (auto& row : elements)
      row.resize(Ncols);
  }

  void resize(size_t Nrows, size_t Ncols, const NumberFormat& val)
  {
    elements.resize(Nrows);
    for (auto& row : elements)
      row.resize(Ncols,val);
  }

  void reserve(size_t Nrows)
  {
    elements.reserve(Nrows);
  }

  void push_back(const std::vector<NumberFormat>& val)
    {elements.push_back(val);}
  void pop_back() {elements.pop_back();}

  bool empty() const noexcept {return elements.empty();}

  //============================================= Iterator access
  typename std::vector<std::vector<NumberFormat>>::iterator
    begin() {return elements.begin();}

  typename std::vector<std::vector<NumberFormat>>::iterator
    end() {return elements.end();}

  size_t size() const {return elements.size();}

  MatDim Dimensions() const
  {
    if (elements.empty())
      return {0,0};
    else
      return {elements.size(),elements[0].size()};
  }

  void bounds_check_rows_cols(const MatDim a, const MatDim b) const
  {
    if ((a.first != b.first) or (a.second != b.second))
    {
      std::length_error excp("Mismatched square sizes of DMatrixNX");
      throw excp;
    }
  }

  void bounds_check_colsA_rowsB(const MatDim a, const MatDim b) const
  {
    if (a.first != b.second)
    {
      std::length_error excp("Mismatched matrix A rows with matrix B cols"
                             " in DMatrixNX");
      throw excp;
    }
  }

  //============================================= Addition
  /**Component-wise addition of two matrices.
   * \f$ \vec{w} = \vec{x} + \vec{y} \f$*/
  DMatrixNX operator+(const DMatrixNX& rhs) const
  {
    auto dim = Dimensions();
    bounds_check_rows_cols(dim,rhs.Dimensions());
    DMatrixNX<NumberFormat> newVector(dim.first,dim.second,0.0);
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        newVector.elements[i][j] = elements[i][j] + rhs.elements[i][j];

    return newVector;
  }

  /**In-place component-wise addition of two vectors.
 * \f$ \vec{x} = \vec{x} + \vec{y} \f$*/
  DMatrixNX& operator+=(const DMatrixNX& rhs)
  {
    auto dim = Dimensions();
    bounds_check_rows_cols(dim,rhs.Dimensions());
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        elements[i][j] = elements[i][j] + rhs.elements[i][j];

    return *this;
  }

  //=========================================== Subtraction
  /**Component-wise subtraction.
 * \f$ \vec{w} = \vec{x} - \vec{y} \f$*/
  DMatrixNX operator-(const DMatrixNX& rhs) const
  {
    auto dim = Dimensions();
    bounds_check_rows_cols(dim,rhs.Dimensions());
    DMatrixNX<NumberFormat> newVector(dim.first,dim.second,0.0);
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        newVector.elements[i][j] = elements[i][j] - rhs.elements[i][j];

    return newVector;
  }

  /**In-place component-wise subtraction.
 * \f$ \vec{x} = \vec{x} - \vec{y} \f$*/
  DMatrixNX& operator-=(const DMatrixNX& rhs)
  {
    auto dim = Dimensions();
    bounds_check_rows_cols(dim,rhs.Dimensions());
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        elements[i][j] = elements[i][j] - rhs.elements[i][j];

    return *this;
  }

  //=========================================== Multiplication
  /**Vector component-wise multiplication by scalar.
 * \f$ \vec{w} = \vec{x} \alpha \f$*/
  DMatrixNX operator*(const NumberFormat value) const
  {
    auto dim = Dimensions();
    DMatrixNX<NumberFormat> newVector(dim.first,dim.second);
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        newVector.elements[i][j] = elements[i][j]*value;

    return newVector;
  }

  /**Vector in-place component-wise multiplication by scalar.
 * \f$ \vec{x} = \vec{x} \alpha \f$*/
  DMatrixNX& operator*=(const NumberFormat value)
  {
    auto dim = Dimensions();
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        elements[i][j] *= value;

    return *this;
  }

  /** Matrix-Matrix multiplication */
  DMatrixNX operator*(const DMatrixNX& rhs)
  {
    auto dimA = Dimensions();
    auto dimB = rhs.Dimensions();
    bounds_check_colsA_rowsB(dimA,dimB);

    DMatrixNX<NumberFormat> newMatrix(dimA.first,dimB.second);
    unsigned int istar=0;
    unsigned int jstar=0;
    for (unsigned int i=0; i<dimA.first; ++i)
    {
      for (unsigned int j=0; j<dimA.first; ++j)
      {
        NumberFormat value = 0.0;
        for (unsigned int k=0; k<dimA.second; ++k)
          value += elements[i][k]*rhs.elements[k][j];

        newMatrix[istar][jstar] = value;
        ++jstar;
      }
      ++istar;
      jstar = 0;
    }//for i

    return newMatrix;
  }

  /** Matrix-Vector multiplication */
  DVectorNX<NumberFormat> operator*(DVectorNX<NumberFormat>& V)
  {
    auto dimA = Dimensions();
    auto dimV = V.size();

    if (dimA.second != dimV)
    {
      std::length_error excp("Mismatched matrix vector sizes in"
                             " matrix-vector multiplication: DMatrixNX");
      throw excp;
    }

    DVectorNX<NumberFormat> newV(dimA.first);
    unsigned int k=0;
    for (unsigned int i=0; i<dimA.first; ++i)
    {
      NumberFormat value = 0.0;
      for (unsigned int j=0; j<dimA.second; ++j)
        value += elements[i][j]*V[j];
      newV[k] = value;
      ++k;
    }

    return newV;
  }

  //=========================================== Division
  /**Vector component-wise division by scalar.
 * \f$ w_i = \frac{x_i}{\alpha} \f$*/
  DMatrixNX operator/(const NumberFormat value) const
  {
    auto dim = Dimensions();
    DMatrixNX<NumberFormat> newVector(dim.first,dim.second);
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        newVector.elements[i][j] = elements[i][j]/value;

    return newVector;
  }

  /**Vector in-place component-wise division by scalar.
 * \f$ x_i = \frac{x_i}{\alpha} \f$*/
  DMatrixNX& operator/=(const NumberFormat value)
  {
    auto dim = Dimensions();
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        elements[i][j] /= value;

    return *this;
  }

  //============================================= Operations
  /** Set the diagonal using a vector.*/
  void SetDiagonal(DVectorNX<NumberFormat>& V)
  {
    auto dimA = Dimensions();
    auto dimV = V.size();

    if ((dimA.first != dimV) or (dimA.second != dimV))
    {
      std::length_error excp("Mismatched matrix vector sizes in"
                             " matrix-vector diagonal assignment: DMatrixNX");
      throw excp;
    }

    for (int i=0; i<dimA.first; ++i)
      elements[i][i] = V[i];
  }

  /** Set the diagonal using value.*/
  void SetDiagonal(NumberFormat val)
  {
    auto dimA = Dimensions();

    for (int i=0; i<dimA.first; ++i)
      elements[i][i] = val;
  }
};

/**Multiplication by a scalar from the left.*/
template<class NumberFormat>
chi_math::DMatrixNX<NumberFormat>
operator*(const double value, chi_math::DMatrixNX<NumberFormat>& that)
{
  auto dim = that.Dimensions();
  chi_math::DMatrixNX<NumberFormat> newMatrix(dim.first,dim.second);
  for (int i=0; i<dim.first; ++i)
    for (int j=0; j<dim.second;++j)
      newMatrix.elements[i][j] = that.elements[i][j]*value;

  return newMatrix;
}

#endif