#ifndef CHI_MATH_LINEAR_INDEXER_H
#define CHI_MATH_LINEAR_INDEXER_H

#include "../chi_math_incdef.h"
#include "../CHI_MATH_MATRIX/chi_math_matrix.h"
#include "../../CHI_VECTOR/chi_vector.h"

class CHI_MATH_LINEAR_INDEXER;
//template class CHI_VECTOR < CHI_MATH_LINEAR_INDEXER >;

//######################################################### Sub-class Definition
class CHI_MATH_LINEAR_INDEXER
{
public:
	long int minIndex=0;
	long int maxIndex=0;
	double   minValue=0;
	double   maxValue=0;
	int      level=-1;
	bool     nextLevelDefined=false;
	CHI_VECTOR<CHI_MATH_LINEAR_INDEXER>     indexers;

public:
	CHI_MATH_LINEAR_INDEXER();
	void     CreateRecursion(CHI_MATH_MATRIX* matrix, int levelAmount, int intervalsPerLevel);
	void     FindIndexRanges(double lookUpValue, long int& kmin, long int& kmax,CHI_MATH_MATRIX* matrix);
};

#endif
