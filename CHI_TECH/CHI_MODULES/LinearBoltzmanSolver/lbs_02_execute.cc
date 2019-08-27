#include "lbs_linear_boltzman_solver.h"
#include "CHI_MODULES/CHI_NPTRANSPORT/IterativeMethods/lbs_iterativemethods.h"

#include <chi_mpi.h>
#include <chi_log.h>
#include <CHI_CONSOLE/chi_console.h>

extern CHI_MPI     chi_mpi;
extern CHI_LOG     chi_log;
extern CHI_CONSOLE chi_console;

//###################################################################
/**Execute the solver.*/
void LinearBoltzmanSolver::Execute()
{
  MPI_Barrier(MPI_COMM_WORLD);
  for (int gs=0; gs<group_sets.size(); gs++)
  {
    chi_log.Log(LOG_0)
      << "\n********* Initializing Groupset " << gs << "\n" << std::endl;

    group_sets[gs]->BuildDiscMomOperator(options.scattering_order);
    group_sets[gs]->BuildMomDiscOperator(options.scattering_order);
    group_sets[gs]->BuildSubsets();

    ComputeSweepOrderings(group_sets[gs]);
    InitFluxDataStructures(group_sets[gs]);

    InitWGDSA(group_sets[gs]);
    InitTGDSA(group_sets[gs]);

    SolveGroupset(gs);

    ResetSweepOrderings(group_sets[gs]);

    MPI_Barrier(MPI_COMM_WORLD);
  }



  chi_log.Log(LOG_0) << "NPTransport solver execution completed\n";
}


//###################################################################
/**Solves a single groupset.*/
void LinearBoltzmanSolver::SolveGroupset(int group_set_num)
{
  NPT_GROUPSET* group_set = group_sets[group_set_num];
  if (group_set->iterative_method == NPT_CLASSICRICHARDSON)
  {
    ClassicRichardson(group_set_num);
  }
  else if (group_set->iterative_method == NPT_GMRES)
  {
    GMRES(group_set_num);
  }

  chi_log.Log(LOG_0)
    << "Groupset solve complete.                  Process memory = "
    << std::setprecision(3)
    << chi_console.GetMemoryUsageInMB() << " MB";
}

