#ifndef CHI_HYDROVOLUME_H
#define CHI_HYDROVOLUME_H

#include "../hydro_component.h"
//######################################################### CLASS DEF
/**Object for controlling hydrodynamic volumes.

\author Jan*/
class HYDRO_VOLUME: public HYDRO_COMPONENT
{
public:
	//int             componentType;  ///< Type code
	double          L;              ///< Length [m]
	double          theta;          ///< Cartesian angle [deg]
	double          phi;            ///< Azimuthal angle [deg]
	double          A;              ///< Cross-sectional area [m^2]
	double          Dh;             ///< Hydraulic diameter [m]
	double          roughness;      ///< Wall roughness [m]
	//int             index;
	
	double          location[3];    ///< Location of initial point
	double          endpoint[3];    ///< End-point
	
	int             leftJunctionIndex;
	int             rigtJunctionIndex;
	
	double          P;              ///< Pressure [Pa]
	double          T_g;            ///< Vapor Temperature [K]
	double          T_f;            ///< Fluid Temperature [K]
	double          rho_g;          ///< Vapor Density [kg/m^3]
	double          rho_f;          ///< Fluid Density [kg/m^3]
	double          mu_g;           ///< Vapor Viscosity [kg/ms]
	double          mu_f;           ///< Fluid Viscosity [kg/ms]
	double          alpha_g;        ///< Vapor volume fraction
	double          alpha_f;        ///< Vapor volume fraction
	double          v_g;            ///< Vapor Velocity [m/s]
	double          v_f;            ///< Fluid Velocity [m/s]
	
public:
					HYDRO_VOLUME();
};

#endif