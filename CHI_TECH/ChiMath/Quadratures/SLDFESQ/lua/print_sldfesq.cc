#include "ChiLua/chi_lua.h"

#include "../sldfe_sq.h"

#include "ChiMath/chi_math.h"
extern ChiMath&     chi_math_handler;

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/** Outputs the quadrature information to python format.
\param handle int Handle to the reference quadrature.
\param file_name_prefix string Prefix to be used in front of file.

##_

###Example:
Example of printing a quadrature:
Example with refinement level 2 and a triple directional refinement:
\code
pquad = chiCreateSLDFESQAngularQuadrature(2)
chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},45.0*math.pi/180,false)
chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},23.0*math.pi/180,false)
chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},12.0*math.pi/180,false)
chiPrintToPythonSLDFESQAngularQuadrature(pquad,"YQuad_");
\endcode

\image html "SLDFESQr.png" width=500px

\ingroup LuaSLDFESQ
\author Jan */
int chiPrintToPythonSLDFESQAngularQuadrature(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("chiPrintToPythonSLDFESQAngularQuadrature",2,num_args);

  int handle = lua_tonumber(L,1);
  const char* file_name = lua_tostring(L,2);

  try{
    auto ref_quadrature = chi_math_handler.angular_quadratures.at(handle);
    if (ref_quadrature->type == chi_math::AngularQuadratureType::SLDFESQ)
    {
      auto sldfesq = std::dynamic_pointer_cast<
        chi_math::SimplifiedLDFESQ::Quadrature>(ref_quadrature);

      if (chi_mpi.location_id == 0)
      {
        sldfesq->output_filename_prefix = file_name;
        sldfesq->PrintQuadratureToFile();
      }
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "chiPrintToPythonSLDFESQAngularQuadrature: "
           "Invalid angular quadrature type.";
      exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiPrintToPythonSLDFESQAngularQuadrature: "
         "Invalid handle to angular quadrature.";
    exit(EXIT_FAILURE);
  }
  catch (...)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiPrintToPythonSLDFESQAngularQuadrature: "
         "Call failed with unknown error.";
    exit(EXIT_FAILURE);
  }


  return 0;
}