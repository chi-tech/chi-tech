#include "lbs_linear_boltzman_solver.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog chi_log;
extern ChiMPI chi_mpi;

//###################################################################
/**Prints header information of simulation.*/
void LinearBoltzman::Solver::PrintSimHeader()
{
  if (chi_mpi.location_id == 0)
  {
    chi_log.Log(LOG_0) << "Initializing LBS Solver\n";
    chi_log.Log(LOG_0) << "Scattering order    : "
                       << options.scattering_order << std::endl;
    chi_log.Log(LOG_0) << "Number of Groups    : "
                       << groups.size() << std::endl;
    chi_log.Log(LOG_0) << "Number of Group sets: "
                       << group_sets.size() << std::endl;
    chi_log.Log(LOG_0) << "Number of Moments   : "
                       << num_moments << std::endl;

    //================================================== Output Groupsets
    for (int gs=0; gs<group_sets.size(); gs++)
    {
      char buf_pol[20];
      std::string outstr;
      int counter = 0;

      chi_log.Log(LOG_0) << "\n***** Groupset " << gs << " *****\n";
      chi_log.Log(LOG_0) << "Groups: ";
      outstr = std::string("");
      counter = 0;
      for (auto group : group_sets[gs]->groups)
      {
        snprintf(buf_pol,20,"%5d ",group->id);
        outstr += std::string(buf_pol);
        counter++;
        if (counter == 12)
        {
          counter = 0;
          chi_log.Log(LOG_0) << outstr << "\n";
          outstr = std::string("");
        }

      }//for g
      chi_log.Log(LOG_0) << outstr << "\n\n";


      chi_math::ProductQuadrature* quad = group_sets[gs]->quadrature;
      chi_log.Log(LOG_0VERBOSE_1) << "Product Quadrature polar angles:\n";
      outstr = std::string("");
      counter = 0;
      for (auto polar_angle : quad->polar_ang)
      {
        snprintf(buf_pol,20,"%9.3f ",polar_angle);
        outstr += std::string(buf_pol);
        counter++;
        if (counter == 6)
        {
          counter = 0;
          chi_log.Log(LOG_0VERBOSE_1) << outstr << "\n";
          outstr = std::string("");
        }
      }
      chi_log.Log(LOG_0VERBOSE_1) << outstr << "\n\n";


      chi_log.Log(LOG_0VERBOSE_1) << "Product Quadrature azimuthal angles:\n";
      outstr = std::string("");
      counter = 0;
      for (auto azimu_ang : quad->azimu_ang)
      {
        snprintf(buf_pol,20,"%9.3f ",azimu_ang);
        outstr += std::string(buf_pol);
        counter++;
        if (counter == 6)
        {
          counter = 0;
          chi_log.Log(LOG_0VERBOSE_1) << outstr << "\n";
          outstr = std::string("");
        }
      }
      chi_log.Log(LOG_0VERBOSE_1) << outstr << "\n\n";
    }
  }
}