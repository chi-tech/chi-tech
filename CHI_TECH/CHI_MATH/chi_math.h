#ifndef CHI_MATH_H
#define CHI_MATH_H

#include "chi_math_incdef.h"

#include "SparseMatrix/chi_math_sparse_matrix_v1.h"

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
public:
	//00 Constructor
						CHI_MATH();
	//01 Utility

	//02 Linear system solvers
	void                GaussElimination(std::vector<std::vector<double> >& A, std::vector<double>& b, int n);
};

namespace chi_math
{
  class SparseMatrix;

  class CDFSampler;
  int SampleCDF(double x, std::vector<double> cdf_bin);
}


#endif
