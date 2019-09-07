#ifndef CHI_PHYSICS_H
#define CHI_PHYSICS_H

#include<iostream>


#include "CHI_SOLVER/chi_solver.h"
#include "CHI_PHYSICSMATERIAL/chi_physicsmaterial.h"
#include "CHI_FIELDFUNCTION/chi_fieldfunction.h"


//############################################################################# CLASS DEF
/** Object for controlling real-time physics.*/
class CHI_PHYSICS
{
	public:
	double    					physicsTimestep;
	double              performanceData[10000];     ///< Misc performance data;
	double              physicsTimeCost;

  std::vector<chi_physics::Solver*>        solver_stack;
  std::vector<chi_physics::Material*>      material_stack;
  std::vector<chi_physics::FieldFunction*> fieldfunc_stack;

	public:
	//00
			CHI_PHYSICS();
	int  InitPetSc(int argc, char** argv);
	//01
	void	RunPhysicsLoop();
	//02
	void    PrintPerformanceData(char* fileName);
};




#endif
