#ifndef CHI_PHYSICS_H
#define CHI_PHYSICS_H

#include<iostream>


#include "SolverBase/chi_solver.h"
#include "PhysicsMaterial/chi_physicsmaterial.h"
#include "ChiPhysics/FieldFunction/fieldfunction.h"


//############################################################################# CLASS DEF
/** Object for controlling real-time physics.*/
class ChiPhysics
{
public:
	double    					physicsTimestep;
	double              performanceData[10000];     ///< Misc performance data;
	double              physicsTimeCost;

  std::vector<chi_physics::Solver*>        solver_stack;
  std::vector<chi_physics::Material*>      material_stack;
  std::vector<chi_physics::TransportCrossSections*> trnsprt_xs_stack;
  std::vector<chi_physics::FieldFunction*> fieldfunc_stack;

private:
  static ChiPhysics instance;

private:
	//00
			ChiPhysics() noexcept;
public:
  static ChiPhysics& GetInstance()
  {return instance;}
	int  InitPetSc(int argc, char** argv);
	//01
	void	RunPhysicsLoop();
	//02
	void    PrintPerformanceData(char* fileName);

};




#endif
