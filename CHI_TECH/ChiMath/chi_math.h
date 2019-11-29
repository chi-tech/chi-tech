#ifndef CHI_MATH_H
#define CHI_MATH_H

#include "chi_math_incdef.h"


#include "Quadratures/quadrature.h"
#include "Quadratures/product_quadrature.h"


/**\defgroup LuaMath B Math*/

typedef std::vector<double> VecDbl;
typedef std::vector<VecDbl> MatDbl;

//######################################################### Class definition
/** This object handles the flow of data for math.

More information of solvers can be obtained from:
[HERE](http://www.maa.org/press/periodicals/loci/joma/iterative-methods-for-solving-iaxi-ibi-the-sor-method)
*/
class ChiMath
{
public:
	std::vector<chi_math::Quadrature*> quadratures;
	std::vector<chi_math::ProductQuadrature*> product_quadratures;
public:
	//00 Constructor
	ChiMath() noexcept;
	//01 Utility

};

namespace chi_math
{
  class SparseMatrix;

  class CDFSampler;
  int SampleCDF(double x, std::vector<double> cdf_bin);

  //02 Vector operations
  void   PrintVector(const VecDbl& x);
  void   Scale(VecDbl& x, const double& val);
  VecDbl VecMul(const VecDbl& x, const double& val);
  double Vec1Norm(const VecDbl& x);
  double Vec2Norm(const VecDbl& x);
  double VecInfinityNorm(const VecDbl& x);
  double VecPNorm(const VecDbl& x, const double& p);
  double Dot(const VecDbl& x, const VecDbl& y);

  //03 Matrix operations
  void   PrintMatrix(const MatDbl& A);
  void   Scale(MatDbl& A, const double& val);
  MatDbl Transpose(const MatDbl& A);
  void   SwapRow(size_t r1, size_t r2, MatDbl& A);
  void   SwapColumn(size_t c1, size_t c2, MatDbl& A);
  MatDbl MatMul(const MatDbl& A, const double c);
  VecDbl MatMul(const MatDbl& A, const VecDbl& x);
  MatDbl MatMul(const MatDbl& A, const MatDbl& B);
  MatDbl MatAdd(const MatDbl& A, const MatDbl& B);
  MatDbl MatSubtract(const MatDbl& A, const MatDbl& B);
  double Determinant(const MatDbl& A);
  MatDbl SubMatrix( const size_t r,
                    const size_t c,
                    const MatDbl& A );
  void   GaussElimination(MatDbl& A, VecDbl& b, int n);
  MatDbl InverseGEPivoting(const MatDbl& A);
  MatDbl Inverse(const MatDbl& A);

  double PowerIteration(const MatDbl& A,
                        VecDbl& e_vec, int max_it = 2000, double tol = 1.0e-13);
}


#endif
