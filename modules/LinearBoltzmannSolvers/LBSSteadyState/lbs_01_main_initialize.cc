#include "lbs_linear_boltzmann_solver.h"

#include "ChiMPI/chi_mpi.h"
#include "chi_log.h"

//###################################################################
/** Initialize the solver.*/
void lbs::SteadyStateSolver::Initialize()
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

  // Initialize source func
  using namespace std::placeholders;
  active_set_source_function_ =
    std::bind(&SteadyStateSolver::SetSource, this, _1, _2, _3, _4);

  //================================================== Initialize groupsets_ for
  //                                                   sweeping
  chi::log.Log() << "Initializing groupset sweeping data" << TextName() << ".";
  for (auto& groupset : groupsets_)
  {
    ComputeSweepOrderings(groupset);
    InitFluxDataStructures(groupset);

    InitWGDSA(groupset);
    InitTGDSA(groupset);
  }

  InitializeSolverSchemes();           //j
  source_event_tag_ = chi::log.GetRepeatingEventTag("Set Source");
}
