#include "chi_math.h"
#include <assert.h>

//############################################################################# GAUSS ELIMINATION
/** Gauss Elimination without pivoting.*/
void CHI_MATH::GaussElimination(std::vector<std::vector<double> > &A,
	                              std::vector<double> &b, int n)
{
	// Forward elimination
	for(int i = 0;i < n-1;++i)
	{
		const std::vector<double>& ai = A[i];
		double bi = b[i];
		double factor = 1.0/A[i][i];
		for(int j = i+1;j < n;++j)
		{
			std::vector<double>& aj = A[j];
			double val = aj[i] * factor;
			b[j] -= val * bi;
			for(int k = i+1;k < n;++k)
				aj[k] -= val * ai[k];
		}
	}

	// Back substitution
	for(int i = n-1;i >= 0;--i)
	{
		const std::vector<double>& ai = A[i];
		double bi = b[i];
		for(int j = i+1;j < n;++j)
			bi -= ai[j] * b[j];
		b[i] = bi/ai[i];
	}
}














