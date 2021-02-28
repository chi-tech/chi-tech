#ifndef CHI_THERMOALPHA_H
#define CHI_THERMOALPHA_H

#define TA_PRESSURE               1
#define TA_BCPRESSURE             2
#define TA_TEMPERATURE            3
#define TA_BCTEMPERATURE          4
#define TA_VOIDFRACTION           5
#define TA_BCVOIDFRACTION         6
#define TA_VELOCITY               7
#define TA_JVELOCITY              8
#define TA_PRESSURETEMPERATURE    9
#define TA_BCPRESSURETEMPERATURE  10
#define TA_JLOSSCOEFFICIENTS      11

#define TA_FLUID_WATER            1001

#include "../../CHI_VECTOR/chi_vector.h"


#include "hydro_component.h"
#include "hydro_sjunction/hydro_sjunction.h"
#include"../../CHI_GRAPHICS/CHI_MATERIAL/CHI_FLOWPROPERTIES/chi_flowproperties.h"





//######################################################### CLASS DEF
/**Object for controlling thermodynamic system.

\author Jan*/
class ThermoAlpha
{
public:
	int fluidOption;
	CHI_VECTOR<HYDRO_COMPONENT>    hydroComponentStack;
	CHI_VECTOR<HYDRO_SJUNCTION> hydroSJunctionStack;
	CHI_FLOWPROPERTIES*         flProps;
	
	CHI_VECTOR<HYDRO_COMPONENT> nodalization;
	
public:
	//00
			ThermoAlpha();
	//01
	bool    Initialize();
	//02
	bool    Step(double dt);
};



#endif