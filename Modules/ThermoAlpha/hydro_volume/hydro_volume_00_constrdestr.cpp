#include "hydro_volume.h"

//######################################################### Default constructor
/** Default constructor.
\author Jan*/
HYDRO_VOLUME::HYDRO_VOLUME()
{
	componentType = 1;
	L           = 1.0;
	theta       = 0.0;
	phi         = 0.0;
	A           = 3.14*0.1*0.1/4.0;  //100mm diameter pipe
	Dh          = 0.1;
	roughness   = 10.0e-6;
	
	location[0] = 0.0;
	location[1] = 0.0;
	location[1] = 0.5;
	
	P           = 100.0e3;
	T_f         = 293.15;
	T_g         = 293.15;
	alpha_f     = 1.0;
	alpha_g     = 0.0;
	v_f         = 0.0;
	v_g         = 0.0;
}