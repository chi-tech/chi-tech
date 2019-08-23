#ifndef CHI_MATH_H
#define CHI_MATH_H

#include "chi_math_incdef.h"

#include "CHI_MATH_MATRIX/chi_math_matrix.h"
#include "SparseMatrix/chi_math_sparse_matrix_v1.h"
#include "CHI_MATH_LINEAR_INTERPOLATION/chi_math_linear_interpolation.h"
#include "CHI_MATH_SPARSE_SOLVER_CG/chi_math_sparse_solver_cg.h"

#include "CHI_MATH_SOLVER_CG/chi_math_solver_cg.h"
#include "Quadratures/quadrature.h"
#include "Quadratures/product_quadrature.h"


#include "../CHI_VECTOR/chi_vector.h"
#include "chi_math_structs.h"

/**\defgroup LuaMath B Math*/

//######################################################### Class definition
/** This object handles the flow of data for math.

More information of solvers can be obtained from:
[HERE](http://www.maa.org/press/periodicals/loci/joma/iterative-methods-for-solving-iaxi-ibi-the-sor-method)
*/
class CHI_MATH
{
public:
	std::vector<CHI_QUADRATURE*> quadratures;
	std::vector<CHI_PRODUCT_QUADRATURE*> product_quadratures;
	CHI_MATH_MATRIX     IVP_coeff1;
	CHI_MATH_MATRIX     IVP_coeff2;
public:
	//00 Constructor
						CHI_MATH();
	//01 Utility
	double				LinearInterpolate(double lookUpValue, CHI_MATH_MATRIX* lookUpMatrix, int lookUpColumn, int indexColumn = 1);
	void				ReIndexTable(CHI_MATH_MATRIX*& A, CHI_MATH_MATRIX* B);
	void                CreateTestMatrixA(CHI_MATH_MATRIX*& A, CHI_MATH_MATRIX*& b, long int size, bool linear = false, double K = 2.0);
	void                CreateTestMatrixB(CHI_MATH_MATRIX*& A, CHI_MATH_MATRIX*& b, long int size);
	void                CreateTestMatrixC(CHI_MATH_MATRIX*& A, CHI_MATH_MATRIX*& b, long int size);
	//02 Linear system solvers
	void                GaussElimination(std::vector<std::vector<double> >& A, std::vector<double>& b, int n);
	void                TDMA(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* x, CHI_MATH_MATRIX* b, int N);
	long int            SolveUsingJacobiIteration(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b, double tolerance = 1.0e-12, bool initialZero = true);
	long int            SolveUsingSORIteration(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b, double omega = 1.0, double tolerance = 1.0e-12, bool initialZero = true);
	long int            SolveUsingSymSORIteration(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b, double omega = 1.3, double tolerance = 1.0e-12, bool initialZero = true);
	long int            SolveUsingSteepestDescent(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b, double tolerance = 1.0e-12, bool initialZero = true);
	long int            SolveUsingConjugantGrad(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b, double tolerance = 1.0e-12, bool initialZero = true);
	long int            SolveSparseUsingConjugantGrad(CHI_MATH_SPARSE_MATRIX* A, CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b, double tolerance = 1.0e-12, bool initialZero = true);

};

namespace chi_math
{
  class SparseMatrix;

  class CDFSampler;
  int SampleCDF(double x, std::vector<double> cdf_bin);
}


#endif
