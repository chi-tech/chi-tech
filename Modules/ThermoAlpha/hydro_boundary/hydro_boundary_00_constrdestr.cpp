#include"hydro_boundary.h"

//######################################################### Default constructor
/** Default constructor

\author Jan*/
HYDRO_BOUNDARY::HYDRO_BOUNDARY()
{
	componentType = 2;
	
	P           = 100.0e3;
	T_f         = 293.15;
	T_g         = 293.15;
	alpha_f     = 1.0;
	alpha_g     = 0.0;
	
}