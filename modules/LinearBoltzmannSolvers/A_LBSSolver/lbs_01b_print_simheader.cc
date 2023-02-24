#include "lbs_solver.h"

#include "ChiLog/chi_log.h"
#include "ChiMPI/chi_mpi.h"

//###################################################################
/**Prints header information of simulation.*/
void lbs::LBSSolver::PrintSimHeader()
{
  if (chi::mpi.location_id == 0)
  {
    std::stringstream outstr;
    outstr << "\nInitializing LBS SteadyStateSolver with name: "
           << TextName() << "\n\n"
           << "Scattering order    : "
           << options_.scattering_order << "\n"
           << "Number of Groups    : "
           << groups_.size() << "\n"
           << "Number of Group sets: "
           << groupsets_.size() << std::endl;

    //================================================== Output Groupsets
    for (const auto& groupset : groupsets_)
    {
      char buf_pol[20];

      outstr << "\n***** Groupset " << groupset.id_ << " *****\n" << "Groups: ";
      int counter = 0;
      for (auto group : groupset.groups)
      {
        snprintf(buf_pol,20,"%5d ",group.id_);
        outstr << std::string(buf_pol);
        counter++;
        if (counter == 12)
        {
          counter = 0;
          outstr << "\n";
        }

      }//for g
      chi::log.Log() << outstr.str() << "\n" << std::endl;
    }//for gs
  }
}