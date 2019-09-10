#ifndef _chi_meshmatrix3x3_h
#define _chi_meshmatrix3x3_h

#include <cmath>
#include <sstream>
#include "chi_meshvector.h"

struct chi_mesh::Matrix3x3
{
  double vals[9];

  /**Creates a zeros matrix of size 3x3.*/
  Matrix3x3()
  {
    vals[0] = 0.0; vals[1] = 0.0; vals[2] = 0.0;
    vals[3] = 0.0; vals[4] = 0.0; vals[5] = 0.0;
    vals[6] = 0.0; vals[7] = 0.0; vals[8] = 0.0;
  }

  /**Copy constructor*/
  void operator=(const Matrix3x3& inM)
  {
    for (int k=0;k<9;k++)
      this->vals[k] = inM.vals[k];
  }

  /**Matrix addition operator.*/
  Matrix3x3 operator+(const Matrix3x3& inM)
  {
    Matrix3x3 oM;
    for (int k=0;k<9;k++)
      oM.vals[k] = this->vals[k] + inM.vals[k];
    return oM;
  }

  /**Matrix subtraction operator.*/
  Matrix3x3 operator-(const Matrix3x3& inM)
  {
    Matrix3x3 oM;
    for (int k=0;k<9;k++)
      oM.vals[k] = this->vals[k] - inM.vals[k];
    return oM;
  }

  /**Matrix multiply with scalar.*/
  Matrix3x3 operator*(const double value)
  {
    Matrix3x3 oM;
    for (int k=0; k<9; k++)
    {
      oM.vals[k] = this->vals[k]*value;
    }
    return oM;
  }

  /**Matrix multiply with vector.*/
  Vector operator*(const Vector& vec)
  {
    double i_vec[] = {vec.x,vec.y,vec.z};
    double o_vec[] = {0.0,0.0,0.0};

    for (int i=0;i<3;i++)
    {
      for (int j=0;j<3;j++)
      {
        o_vec[i] += this->GetIJ(i,j)*i_vec[j];
      }
    }
    Vector oV(o_vec[0],o_vec[1],o_vec[2]);
    return oV;
  }

  /**Set value at row i and column j.*/
  void SetIJ(int i, int j, double value)
  {
    int k = j + 3*i;
    vals[k] = value;
  }

  /**Add value to value at row i and column j.*/
  void AddIJ(int i, int j, double value)
  {
    int k = j + 3*i;
    vals[k] += value;
  }

  /**Obtain a copy of the value at row i and column j.*/
  double GetIJ(int i, int j)
  {
    int k = j + 3*i;
    return vals[k];
  }

  /**Set row i using a vector.*/
  void SetRowIVec(int i, Vector vec)
  {
    vals[0 + 3*i] = vec.x;
    vals[1 + 3*i] = vec.y;
    vals[2 + 3*i] = vec.z;
  }

  /**Set column j using a vector.*/
  void SetColJVec(int j, Vector vec)
  {
    vals[j + 3*0] = vec.x;
    vals[j + 3*1] = vec.y;
    vals[j + 3*2] = vec.z;
  }

  /**Sets the diagonal of the matrix.*/
  void SetDiagonalVec(double a00,double a11, double a22)
  {
    vals[0 + 3*0] = a00;
    vals[1 + 3*1] = a11;
    vals[2 + 3*2] = a22;
  }

  /**Get the determinant using specified row [default:0].*/
  double Det(int row=0)
  {
    double det=0.0;

    int sign = -1;
    for (int j=0; j<3; j++)
    {
      int k = j + 3*row;
      sign *= -1;

      det += sign*GetIJ(row,j)*MinorIJ(row,j);
    }

    return det;
  }

  /**Get the minor value associated with row ir and column jr.*/
  double MinorIJ(int ir, int jr)
  {
    double a[] = {0.0,0.0,0.0,0.0};

    int k = -1;
    for (int i=0; i<3; i++)
    {
      if (i==ir) continue;

      for (int j=0; j<3; j++)
      {
        if (j==jr) continue;

        k++;
        int kr = j + 3*i;
        a[k] = vals[kr];
      }
    }

    return a[0]*a[3]-a[2]*a[1];
  }

  /**Compute the matrix transpose.*/
  Matrix3x3 Transpose()
  {
    Matrix3x3 oM;
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<3; j++)
      {
        int kO = j + 3*i;
        int kI = i + 3*j;

        oM.vals[kO] = vals[kI];
      }
    }
    return oM;
  }

  /**Compute the matrix transpose.*/
  Matrix3x3 Inverse()
  {
    Matrix3x3 oM;
    Matrix3x3 oMT;

    //================================= Compute matrix of minors
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<3; j++)
      {
        oM.SetIJ(i,j,MinorIJ(i,j));
      }
    }

    //================================= Compute matrix of cofactors
    int sign = -1;
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<3; j++)
      {
        sign*=-1;
        oM.SetIJ(i,j,oM.GetIJ(i,j)*sign);
      }
    }

    //================================= Compute the transpose
    oMT = oM.Transpose();

    //================================= Get determinant
    double det = Det();

    return oMT*(1.0/det);
  }

  /**Outputs the matrix to a stream.*/
  friend std::ostream & operator<< (std::ostream& out, Matrix3x3& inM)
  {
    out << "[";
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<3;j++)
      {
        out << inM.GetIJ(i,j) << " ";
      }
      if (i!=2)
        out << "\n";
      else
        out << "]\n";
    }

    return out;
  }

  std::string PrintS()
  {
    std::stringstream out;
    out << "[";
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<3;j++)
      {
        out << GetIJ(i,j) << " ";
      }
      if (i!=2)
        out << "\n ";
      else
        out << "]";
    }

    return out.str();
  }

};

#endif