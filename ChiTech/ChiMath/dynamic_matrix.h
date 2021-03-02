#ifndef CHI_MATH_DYNAMIC_MATRIX_H
#define CHI_MATH_DYNAMIC_MATRIX_H

#include <vector>
#include <stdexcept>
#include <sstream>

namespace chi_math
{
  template<class NumberFormat>
  class DynamicMatrix;
}

#include "dynamic_vector.h"

//###################################################################
/**General dynamic matrix utility.*/
template<class NumberFormat>
class chi_math::DynamicMatrix
{
  typedef std::pair<size_t,size_t> MatDim;
public:
  std::vector<std::vector<NumberFormat>> elements;

  /**Default constructor. Does nothing.*/
  DynamicMatrix()
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are "
                  "supported for DynamicMatrix." );
  }

  /**Constructor with number of entries. Value defaults.*/
  DynamicMatrix(size_t Nrows, size_t Ncols) :
    elements(Nrows,std::vector<NumberFormat>(Ncols))
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are "
                  "supported for DynamicMatrix." );
  }

  /**Constructor with number of entries and default value.*/
  DynamicMatrix(size_t Nrows, size_t Ncols, NumberFormat value) :
    elements(Nrows, std::vector<NumberFormat>(Ncols,value))
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are "
                  "supported for DynamicMatrix." );
  }

  /**Copy constructor.*/
  DynamicMatrix(const DynamicMatrix& other) { elements = other.elements;}

  /**Assignment operator.*/
  DynamicMatrix& operator=(const DynamicMatrix& other)
  {
    elements = other.elements;
    return *this;
  }

  /**Move constructor.*/
  DynamicMatrix(DynamicMatrix&& other) noexcept
  { elements = std::move(other.elements);}

  /**Move assignment operator.*/
  DynamicMatrix& operator=(DynamicMatrix&& other) noexcept
  {
    elements = std::move(other.elements);
    return *this;
  }

  /**Constructor with vector.*/
  explicit
  DynamicMatrix(const std::vector<std::vector<double>>& in) { elements = in;}

  /**Copy constructor with vector.*/
  DynamicMatrix& operator=(const std::vector<std::vector<double>>& in)
  {
    elements = in;
    return *this;
  }

  /**Constructor with vector.*/
  DynamicMatrix(std::initializer_list<std::initializer_list<NumberFormat>> in)
  {
    elements.clear();
    for (auto& v : in)
      elements.push_back(v);
  }

  /**Copy constructor with vector.*/
  DynamicMatrix& operator=(std::initializer_list<std::initializer_list<NumberFormat>> in)
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
    {
      row.resize(Ncols,val);
      for (auto& entry : row)
        entry = val;
    }
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
      throw std::length_error("Mismatched square sizes of DynamicMatrix");
  }

  void bounds_check_colsA_rowsB(const MatDim a, const MatDim b) const
  {
    if (a.first != b.second)
      throw std::length_error("Mismatched matrix A rows with matrix B cols"
                              " in DynamicMatrix");
  }

  //============================================= Addition
  /**Component-wise addition of two matrices.
   * \f$ \vec{w} = \vec{x} + \vec{y} \f$*/
  DynamicMatrix operator+(const DynamicMatrix& rhs) const
  {
    auto dim = Dimensions();
    bounds_check_rows_cols(dim,rhs.Dimensions());
    DynamicMatrix<NumberFormat> newVector(dim.first, dim.second, 0.0);
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        newVector.elements[i][j] = elements[i][j] + rhs.elements[i][j];

    return newVector;
  }

  /**In-place component-wise addition of two vectors.
 * \f$ \vec{x} = \vec{x} + \vec{y} \f$*/
  DynamicMatrix& operator+=(const DynamicMatrix& rhs)
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
  DynamicMatrix operator-(const DynamicMatrix& rhs) const
  {
    auto dim = Dimensions();
    bounds_check_rows_cols(dim,rhs.Dimensions());
    DynamicMatrix<NumberFormat> newVector(dim.first, dim.second, 0.0);
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        newVector.elements[i][j] = elements[i][j] - rhs.elements[i][j];

    return newVector;
  }

  /**In-place component-wise subtraction.
 * \f$ \vec{x} = \vec{x} - \vec{y} \f$*/
  DynamicMatrix& operator-=(const DynamicMatrix& rhs)
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
  DynamicMatrix operator*(const NumberFormat value) const
  {
    auto dim = Dimensions();
    DynamicMatrix<NumberFormat> newVector(dim.first, dim.second);
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        newVector.elements[i][j] = elements[i][j]*value;

    return newVector;
  }

  /**Vector in-place component-wise multiplication by scalar.
 * \f$ \vec{x} = \vec{x} \alpha \f$*/
  DynamicMatrix& operator*=(const NumberFormat value)
  {
    auto dim = Dimensions();
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        elements[i][j] *= value;

    return *this;
  }

  /** Matrix-Matrix multiplication */
  DynamicMatrix operator*(const DynamicMatrix& rhs)
  {
    auto dimA = Dimensions();
    auto dimB = rhs.Dimensions();
    bounds_check_colsA_rowsB(dimA,dimB);

    DynamicMatrix<NumberFormat> newMatrix(dimA.first, dimB.second);
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
  DynamicVector<NumberFormat> operator*(DynamicVector<NumberFormat>& V)
  {
    auto dimA = Dimensions();
    auto dimV = V.size();

    if (dimA.second != dimV)
      throw std::length_error("Mismatched matrix vector sizes in"
                              " matrix-vector multiplication: DynamicMatrix");

    DynamicVector<NumberFormat> newV(dimA.first);
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
  DynamicMatrix operator/(const NumberFormat value) const
  {
    auto dim = Dimensions();
    DynamicMatrix<NumberFormat> newVector(dim.first, dim.second);
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        newVector.elements[i][j] = elements[i][j]/value;

    return newVector;
  }

  /**Vector in-place component-wise division by scalar.
 * \f$ x_i = \frac{x_i}{\alpha} \f$*/
  DynamicMatrix& operator/=(const NumberFormat value)
  {
    auto dim = Dimensions();
    for (int i=0; i<dim.first; ++i)
      for (int j=0; j<dim.second;++j)
        elements[i][j] /= value;

    return *this;
  }

  //============================================= Operations
  /** Set the diagonal using a vector.*/
  void SetDiagonal(DynamicVector<NumberFormat>& V)
  {
    auto dimA = Dimensions();
    auto dimV = V.size();

    if ((dimA.first != dimV) or (dimA.second != dimV))
    {
      throw std::length_error("Mismatched matrix vector sizes in"
                              " matrix-vector diagonal assignment: DynamicMatrix");
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

  /**Prints the matrix to a string and then returns the string.*/
  std::string PrintStr() const
  {
    auto dim = Dimensions();
    std::stringstream out;

    for (int i = 0; i<dim.first ; ++i)
    {
      for (int j=0; j<(dim.second-1); ++j)
        out<<elements[i][j]<<" ";
      out<<elements[i][dim.second-1];

      if (i<(dim.first -1)) out << "\n";
    }

    return out.str();
  }
};

/**Multiplication by a scalar from the left.*/
template<class NumberFormat>
chi_math::DynamicMatrix<NumberFormat>
operator*(const double value, chi_math::DynamicMatrix<NumberFormat>& that)
{
  auto dim = that.Dimensions();
  chi_math::DynamicMatrix<NumberFormat> newMatrix(dim.first, dim.second);
  for (int i=0; i<dim.first; ++i)
    for (int j=0; j<dim.second;++j)
      newMatrix.elements[i][j] = that.elements[i][j]*value;

  return newMatrix;
}

#endif //CHI_MATH_DYNAMIC_MATRIX_H