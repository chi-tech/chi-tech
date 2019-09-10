#include "lbs_linear_boltzman_solver.h"
#include "ChiConsole/chi_console.h"

#include <ChiMesh/SweepUtilities/chi_FLUDS.h>
#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include <ChiMesh/VolumeMesher/Predefined2D/volmesher_predefined2d.h>
#include <ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h>

#include <chi_log.h>
#include <chi_mpi.h>

#include <iomanip>

extern ChiConsole chi_console;
extern ChiLog     chi_log;
extern ChiMPI     chi_mpi;


//###################################################################
/**Initializes fluds data structures.*/
void LinearBoltzmanSolver::InitFluxDataStructures(LBSGroupset *groupset)
{
  chi_mesh::MeshHandler* handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher* mesher = handler->volume_mesher;

  if ((typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D)) or
      (typeid(*mesher) == typeid(chi_mesh::VolumeMesherExtruder)))
  {
    //================================================== Angle Aggregation
    if      (groupset->angleagg_method == NPT_ANGAGG_SINGLE)
    {
      InitAngleAggSingle(groupset);
    }
    else if (groupset->angleagg_method == NPT_ANGAGG_POLAR)
    {
      InitAngleAggPolar(groupset);
    }
  }
  else
  {
    InitAngleAggSingle(groupset);
  }






  //                   FIRST PASS                   //
  //================================================== Loop over angle set groups
  for (int ag=0; ag<groupset->angle_agg->angle_set_groups.size(); ag++)
  {
    int num_angle_sets = groupset->angle_agg->
                         angle_set_groups[ag]->angle_sets.size();

    //====================================== Loop over anglesets
    for (int as=0; as<num_angle_sets; as++)
    {
      chi_mesh::SweepManagement::AngleSet* angle_set =
        groupset->angle_agg->angle_set_groups[ag]->angle_sets[as];

      angle_set->fluds =
        new chi_mesh::SweepManagement::FLUDS(angle_set->GetNumGrps());

      angle_set->fluds->InitializeAlphaElements(angle_set->GetSPDS());
    }
  }//for ag

  chi_log.Log(LOG_0)
    << "Initialized FLUDS Alpha Elements. ";




  //                  SECOND PASS                   //
  //================================================== Loop over angle set groups
  int anglset_counter=-1;
  for (int ag=0; ag<groupset->angle_agg->angle_set_groups.size(); ag++)
  {
    int num_angle_sets = groupset->angle_agg->
                         angle_set_groups[ag]->angle_sets.size();

    //====================================== Loop over anglesets
    for (int as=0; as<num_angle_sets; as++)
    {
      chi_mesh::SweepManagement::AngleSet* angle_set =
        groupset->angle_agg->angle_set_groups[ag]->angle_sets[as];

      anglset_counter++;
      angle_set->fluds->InitializeBetaElements(angle_set->GetSPDS());
    }
  }//for ag

  chi_log.Log(LOG_0)
    << "Initialized FLUDS Beta Elements. "
    << "         Process memory = "
    << std::setprecision(3) << chi_console.GetMemoryUsageInMB()
    << " MB.";





  MPI_Barrier(MPI_COMM_WORLD);
}