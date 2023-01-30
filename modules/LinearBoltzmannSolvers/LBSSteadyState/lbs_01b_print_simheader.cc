#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/Quadratures/angular_product_quadrature.h"

#include "ChiLog/chi_log.h"
#include "ChiMPI/chi_mpi.h"

;


//###################################################################
/**Prints header information of simulation.*/
void lbs::SteadyStateSolver::PrintSimHeader()
{
  if (chi::mpi.location_id == 0)
  {
    chi::log.Log() << "\nInitializing LBS SteadyStateSolver with name: "
                       << TextName() << "\n\n";
    chi::log.Log() << "Scattering order    : "
                       << options.scattering_order << std::endl;
    chi::log.Log() << "Number of Groups    : "
                       << groups.size() << std::endl;
    chi::log.Log() << "Number of Group sets: "
                       << groupsets.size() << std::endl;

    //================================================== Output Groupsets
    for (int gs=0; gs < groupsets.size(); gs++)
    {
      char buf_pol[20];
      std::string outstr;
      int counter = 0;

      chi::log.Log() << "\n***** Groupset " << gs << " *****\n";
      chi::log.Log() << "Groups: ";
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
          chi::log.Log() << outstr << "\n";
          outstr = std::string("");
        }

      }//for g
      chi::log.Log() << outstr << "\n\n";


      auto quad = groupsets[gs].quadrature;

      if (quad->type == chi_math::AngularQuadratureType::ProductQuadrature)
      {
        auto product_quadrature =
          std::static_pointer_cast<chi_math::ProductQuadrature>(quad);


        chi::log.Log0Verbose1() << "Product Quadrature polar angles:\n";
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
            chi::log.Log0Verbose1() << outstr << "\n";
            outstr = std::string("");
          }
        }
        chi::log.Log0Verbose1() << outstr << "\n\n";


        chi::log.Log0Verbose1() << "Product Quadrature azimuthal angles:\n";
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
            chi::log.Log0Verbose1() << outstr << "\n";
            outstr = std::string("");
          }
        }
        chi::log.Log0Verbose1() << outstr << "\n\n";
      }//if product quadrature




    }
  }
}