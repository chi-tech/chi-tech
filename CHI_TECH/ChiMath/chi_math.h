#ifndef CHI_MATH_H
#define CHI_MATH_H

#include "chi_math_incdef.h"


#include "Quadratures/quadrature.h"
#include "Quadratures/product_quadrature.h"


/**\defgroup LuaMath B Math*/

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
	ChiMath();
	//01 Utility

	//02 Linear system solvers
	void GaussElimination(std::vector<std::vector<double> >& A, std::vector<double>& b, int n);
};

namespace chi_math
{
  class SparseMatrix;

  class CDFSampler;
  int SampleCDF(double x, std::vector<double> cdf_bin);
}


#endif
