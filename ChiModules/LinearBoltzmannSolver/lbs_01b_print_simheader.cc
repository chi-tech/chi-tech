#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/Quadratures/product_quadrature.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;


//###################################################################
/**Prints header information of simulation.*/
void lbs::SteadySolver::PrintSimHeader()
{
  if (chi::mpi.location_id == 0)
  {
    chi_log.Log(LOG_0) << "\nInitializing LBS Solver with name: "
                       << TextName() << "\n\n";
    chi_log.Log(LOG_0) << "Scattering order    : "
                       << options.scattering_order << std::endl;
    chi_log.Log(LOG_0) << "Number of Groups    : "
                       << groups.size() << std::endl;
    chi_log.Log(LOG_0) << "Number of Group sets: "
                       << groupsets.size() << std::endl;

    //================================================== Output Groupsets
    for (int gs=0; gs < groupsets.size(); gs++)
    {
      char buf_pol[20];
      std::string outstr;
      int counter = 0;

      chi_log.Log(LOG_0) << "\n***** Groupset " << gs << " *****\n";
      chi_log.Log(LOG_0) << "Groups: ";
      outstr = std::string("");
      counter = 0;
      for (auto group : groupsets[gs].groups)
      {
        snprintf(buf_pol,20,"%5d ",group.id);
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


      auto quad = groupsets[gs].quadrature;

      if (quad->type == chi_math::AngularQuadratureType::ProductQuadrature)
      {
        auto product_quadrature =
          std::static_pointer_cast<chi_math::ProductQuadrature>(quad);


        chi_log.Log(LOG_0VERBOSE_1) << "Product Quadrature polar angles:\n";
        outstr = std::string("");
        counter = 0;
        for (auto polar_angle : product_quadrature->polar_ang)
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
        for (auto azimu_ang : product_quadrature->azimu_ang)
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
      }//if product quadrature




    }
  }
}