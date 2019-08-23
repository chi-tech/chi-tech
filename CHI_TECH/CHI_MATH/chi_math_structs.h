#ifndef CHI_MATH_STRUCTS_H
#define CHI_MATH_STRUCTS_H

#include<iostream>
#include<vector>

//========================================================= Structs
/** Structure for holding the solution of an IVP.*/
struct CHI_MATH_IVP_SOLUTION
{
	CHI_VECTOR<double> t;
	CHI_VECTOR<double> x;
	CHI_VECTOR<double> e;
};



#endif
