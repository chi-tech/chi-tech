#include "chi_runtime.h"

extern bool sim_option_interactive;

//######################################################### Program entry point
/** Program entry point.

\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.

*/
int main(int argc, char** argv)
{
  ChiTechInitialize(argc,argv);

  if (sim_option_interactive)
    ChiTechRunInteractive(argc, argv);
  else
    ChiTechRunBatch(argc, argv);

  ChiTechFinalize();

  return 0;
}