#include "chi_physics.h"
#include <petscksp.h>

//############################################################################# Default constructor
/** Default constructor.*/
CHI_PHYSICS::CHI_PHYSICS()
{
	this->physicsTimestep=16.66667;
	//this->physicsTimestep=1000;
	this->physicsTimeCost=0;

	for (int k=0;k<10000;k++)
	{
		this->performanceData[k]=0.0;
	}
}

/**Initializes PetSc for use by all entities.*/
int CHI_PHYSICS::InitPetSc(int argc, char** argv)
{
	PetscErrorCode ierr;
	PetscMPIInt    size;
	ierr = PetscInitialize(&argc,&argv,(char*)0,NULL);
	if (ierr) return ierr;

	PetscOptionsSetValue(NULL,"-error_output_stderr",NULL);
  PetscOptionsSetValue(NULL,"-no_signal_handler",NULL);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

	return 0;
}
