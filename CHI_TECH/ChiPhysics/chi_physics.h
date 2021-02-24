#ifndef CHI_PHYSICS_H
#define CHI_PHYSICS_H

#include<iostream>


#include "SolverBase/chi_solver.h"
#include "PhysicsMaterial/chi_physicsmaterial.h"
#include "CHI_TECH/ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"


//############################################################################# CLASS DEF
/** Object for controlling real-time physics.*/
class ChiPhysics
{
public:
	double    					physicsTimestep=16.66667;
	double              physicsTimeCost=0.0;

  std::vector<chi_physics::Solver*>                        solver_stack;
  std::vector<std::shared_ptr<chi_physics::Material>>      material_stack;
  std::vector<std::shared_ptr<chi_physics::TransportCrossSections>> trnsprt_xs_stack;
  std::vector<std::shared_ptr<chi_physics::FieldFunction>> fieldfunc_stack;

private:
  static ChiPhysics instance;

private:
	//00
			ChiPhysics() noexcept;
public:
  static ChiPhysics& GetInstance() noexcept
  {return instance;}
	int  InitPetSc(int argc, char** argv);
	//01
	void	RunPhysicsLoop();
	//02
	void    PrintPerformanceData(char* fileName);

};




#endif
