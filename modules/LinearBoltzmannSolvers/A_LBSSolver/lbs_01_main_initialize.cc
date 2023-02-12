#include "lbs_solver.h"

#include "ChiMPI/chi_mpi.h"
#include "chi_log.h"

//###################################################################
/** Initialize the solver.*/
void lbs::LBSSolver::Initialize()
{
  PerformInputChecks();                //a assigns num_groups_ and grid_ptr_
  PrintSimHeader();                    //b

  MPI_Barrier(MPI_COMM_WORLD);

  InitMaterials();                     //c
  InitializeSpatialDiscretization();   //d
  InitializeGroupsets();               //e
  ComputeNumberOfMoments();            //f
  InitializeParrays();                 //g
  InitializeBoundaries();              //h
  InitializePointSources();            //i

  source_event_tag_ = chi::log.GetRepeatingEventTag("Set Source");
}
