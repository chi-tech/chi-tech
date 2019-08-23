#ifndef CHI_TABLE_H
#define CHI_TABLE_H

#include "../../CHI_MATH/CHI_MATH_MATRIX/chi_math_matrix.h"

//############################################################################# Class Def
/**Object for implementing tables.

\author Jan*/
class CHI_TABLE
{
public:
	char    		name[200];
	char*           columnNames;
	char            fileName[300];
	CHI_MATH_MATRIX rawMatrix;
	
	
public:
	//00
				CHI_TABLE();
};

#endif
