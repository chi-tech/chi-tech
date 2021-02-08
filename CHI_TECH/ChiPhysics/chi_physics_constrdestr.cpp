#include "chi_physics.h"
#include <petscksp.h>

extern bool chi_allow_petsc_error_handler;

//############################################################################# Default constructor
/** Default constructor.*/
ChiPhysics::ChiPhysics() noexcept
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
int ChiPhysics::InitPetSc(int argc, char** argv)
{
	PetscErrorCode ierr;
	PetscMPIInt    size;

  PetscOptionsInsertString(NULL,"-error_output_stderr");
  if (not chi_allow_petsc_error_handler)
    PetscOptionsInsertString(NULL,"-no_signal_handler");
  PetscOptionsInsertString(NULL,"-on_error_abort");

	ierr = PetscInitialize(&argc,&argv,(char*)0,NULL);
	if (ierr) return ierr;



	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

	return 0;
}
