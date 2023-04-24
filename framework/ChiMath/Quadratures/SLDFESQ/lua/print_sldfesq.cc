#include "ChiLua/chi_lua.h"

#include "chi_runtime.h"

#include "../sldfe_sq.h"

#include "chi_log.h"

#include "chi_mpi.h"


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
    auto ref_quadrature = chi::angular_quadrature_stack.at(handle);
    if (ref_quadrature->type_ == chi_math::AngularQuadratureType::SLDFESQ)
    {
      auto sldfesq = std::dynamic_pointer_cast<
        chi_math::SimplifiedLDFESQ::Quadrature>(ref_quadrature);

      if (chi::mpi.location_id == 0)
      {
        sldfesq->output_filename_prefix_ = file_name;
        sldfesq->PrintQuadratureToFile();
      }
    }
    else
    {
      chi::log.LogAllError()
        << "chiPrintToPythonSLDFESQAngularQuadrature: "
           "Invalid angular quadrature type.";
     chi::Exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    chi::log.LogAllError()
      << "chiPrintToPythonSLDFESQAngularQuadrature: "
         "Invalid handle to angular quadrature.";
   chi::Exit(EXIT_FAILURE);
  }
  catch (...)
  {
    chi::log.LogAllError()
      << "chiPrintToPythonSLDFESQAngularQuadrature: "
         "Call failed with unknown error.";
   chi::Exit(EXIT_FAILURE);
  }


  return 0;
}