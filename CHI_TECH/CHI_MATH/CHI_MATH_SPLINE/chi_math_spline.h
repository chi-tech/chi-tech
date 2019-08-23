#ifndef CHI_MATH_SPLINE_H
#define CHI_MATH_SPLINE_H

#include "../chi_math_incdef.h"
#include "../../CHI_VECTOR/chi_vector.h"
#include "../CHI_MATH_MATRIX/chi_math_matrix.h"


//############################################################################# Class definition
/** Cubic spline object.

This object has two internal vectors containing the knot interpolation points and function values. 
It also has a dynamically created matrix vector containing all the z values.

The spline can handle 4 types of boundaries:
-# Free ends (second derivatives zero and both ends)
-# Left fixed
-# Right fixed
-# Both ends fixed

For information on the math consult the whitepaper <a href="../whitepages/Math/SplineInterpolation.pdf" target="_blank"><b>"Spline interpolation"</b></a>
.
\author CHIV*/
class CHI_MATH_SPLINE
{
public:
	int                type;
	double             minBoundary;
	double             maxBoundary;
	CHI_VECTOR<double> knotPoints;
	CHI_VECTOR<double> knotValues;

	CHI_MATH_MATRIX*    z;
public:
	//00 Constr
			CHI_MATH_SPLINE();
	//01 General
	bool    InitializeFromFile(char* fileName);
	bool    InitializeFromTable(CHI_MATH_MATRIX* table,int type,double minBound=0.0,double maxBound=0.0);
	void    BuildNaturalSpline();
	void    BuildLeftFixedSpline();
	void    BuildRightFixedSpline();
	void    BuildFixedSpline();
	double  Value(double x);
	double  Derivative(double x);
};

#endif
