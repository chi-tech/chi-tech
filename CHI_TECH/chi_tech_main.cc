#include "chi_runtime.h"

extern bool chi_sim_option_interactive;

//######################################################### Program entry point
/** Program entry point.

\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.

*/
int main(int argc, char** argv)
{
  ChiTech::Initialize(argc, argv);

  int error_code = 0;
  if (chi_sim_option_interactive)
    error_code = ChiTech::RunInteractive(argc, argv);
  else
    error_code = ChiTech::RunBatch(argc, argv);

  ChiTech::Finalize();

  return error_code;
}