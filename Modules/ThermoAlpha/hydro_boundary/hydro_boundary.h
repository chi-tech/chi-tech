#ifndef CHI_HYDROBOUNDARY_H
#define CHI_HYDROBOUNDARY_H

#include "../hydro_component.h"

//############################################################### CLASS DEF
/**Object for controlling hydrodynamic boundary conditions

\author Jan*/
class HYDRO_BOUNDARY : public HYDRO_COMPONENT
{
public:
	//int             componentType;
	double          A;
	double          L;
	double          location[3];    ///< Location of initial point
	//int             index;
	
	int             junctionIndex;
	
	double          P;              ///< Pressure [Pa]
	double          T_g;            ///< Vapor Temperature [K]
	double          T_f;            ///< Fluid Temperature [K]
	double          rho_g;          ///< Vapor Density [kg/m^3]
	double          rho_f;          ///< Fluid Density [kg/m^3]
	double          mu_g;           ///< Vapor Viscosity [kg/ms]
	double          mu_f;           ///< Fluid Viscosity [kg/ms]
	double          alpha_g;        ///< Vapor volume fraction
	double          alpha_f;        ///< Vapor volume fraction
	
	HYDRO_BOUNDARY();
};



#endif