#include "chi_math_spline.h"



//############################################################################# Default constructor
/** Default constructor. Disables threadprotection of the vectors.*/
CHI_MATH_SPLINE::CHI_MATH_SPLINE()
{
	this->type=0;
	this->knotPoints.threadProtectionEnabled=false;
	this->knotValues.threadProtectionEnabled=false;
}
