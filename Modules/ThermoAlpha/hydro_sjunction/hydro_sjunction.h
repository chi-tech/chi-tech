#ifndef HYDRO_SJUNCTION_H
#define HYDRO_SJUNCTION_H

#include "../hydro_component.h"

#define CONNECTED_TO_NOTHING                    0
#define CONNECTED_TO_HYDRO_VOLUME               1
#define CONNECTED_TO_HYDRO_BOUNDARY             2

//######################################################### CLASS DEF
/** Object for controlling Hydrodynamic Single Junctions.
 *
\author  Jan */
class HYDRO_SJUNCTION
{
public:
	HYDRO_COMPONENT* leftHComponent;
	HYDRO_COMPONENT* rigtHComponent;
	
	double floss;       ///< Forward loss coefficient
	double rloss;       ///< Reverse loss coefficient
	
	double v_g;            ///< Vapor Velocity [m/s]
	double v_f;
public:
	HYDRO_SJUNCTION();
};


#endif