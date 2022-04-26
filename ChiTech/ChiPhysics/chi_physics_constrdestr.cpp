#include "chi_physics.h"
#include <petscksp.h>

#include "chi_runtime.h"

//############################################################################# Default constructor
/** Default constructor.*/
ChiPhysics::ChiPhysics() noexcept
{

}

/**Initializes PetSc for use by all entities.*/
int ChiPhysics::InitPetSc(int argc, char** argv)
{
	PetscErrorCode ierr;
	PetscMPIInt    size;

  PetscOptionsInsertString(NULL,"-error_output_stderr");
  if (not chi::run_time::allow_petsc_error_handler)
    PetscOptionsInsertString(NULL,"-no_signal_handler");
//  PetscOptionsInsertString(NULL,"-on_error_abort");
//TODO: Investigate this, causes cfem methods to fail

	ierr = PetscInitialize(&argc,&argv,(char*)0,NULL);
	if (ierr) return ierr;



	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

	return 0;
}
