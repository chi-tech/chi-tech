#ifndef CHI_MATH_LINEAR_INTERPOLATE_H
#define CHI_MATH_LINEAR_INTERPOLATE_H

#include "../chi_math_incdef.h"
#include "../../CHI_VECTOR/chi_vector.h"
#include "chi_math_linear_indexer.h"
#include "../CHI_MATH_MATRIX/chi_math_matrix.h"

//template class CHI_VECTOR < double > ;
//template class CHI_VECTOR < CHI_MATH_LINEAR_INDEXER >;


//############################################################################# Class definition
/** Blanket object for doing repeated linear interpolations.

When the Optimize routine is called, this object creates an array of 
indexer objects. Each indexer object contains pointers to values of the
index column in the main matrix. In other words, let \f$ a_{i,1} = a_{1,1},...,a_{n,1} \f$ be the values of the 
index column in the main optimization matrix, then when the optimize routine is called,
it creates \f$ x_i \biggr |_{a}^{b} = a_{a,1},...,a_{b,1} \f$ where \f$i\f$ ranges over 
the number of intervals per level. Each indexer object has the option to 
procure another level of indexer. Such that the effective amount of search 
operations can be reduced.

## Efficiency
The efficiency is defined in the index value search. For an unoptimized search routine,
the code needs to start from \f$k=1\f$ through \f$k=n-1\f$ and search for the index values
that lay above and below the lookup value. Therefore, each time the function is called,
\f$n\f$ amount of comparison operations are performed at maximum. With this optimization,
The code will first perform comparisons on the first level, for which the default amount of
intervals \f$x=10\f$, and therefore the first level will have \f$n/10\f$ comparisons after which 
the code will run through the remainder n/10 indexes to find the value of interest. So for 
\f$n>>\f$ large, and multiple levels \f$M\f$, each level will have the following amount of 
search operations:

\f[
	n/10 + (n/10)/10 + ... + n/{10}^M = n \sum_{l=1}^{M} 1/{10}^l
\f]
And therefore, there is a limit to having the level approach:
- \f$x=10\f$, 0.111111 reduction
- \f$x=20\f$, 0.0526 reduction
- \f$x=30\f$, 0.0344 reduction
- \f$x=100\f$, 0.0101 reduction
- \f$x=1000\f$, 0.001001 reduction
- \f$x=10^n\f$, \f$1/{10}^n\f$ reduction
*/
class CHI_MATH_LINEAR_INTERPOLATION
{
public:
	int										intervalsPerLevel;        ///< Amount of intervals into which to split each interval
	int										numberOfLevels;           ///< Number of nested index matrices
	CHI_MATH_MATRIX*						mainMatrix;               ///< Main matrix in which optimization will occur
	CHI_VECTOR<CHI_MATH_LINEAR_INDEXER>     indexers;                 ///< List of indexer objects

public:
					 CHI_MATH_LINEAR_INTERPOLATION();
	double           LinearInterpolate(double lookUpValue);
	void             Optimize(CHI_MATH_MATRIX* matrix);
};


#endif
