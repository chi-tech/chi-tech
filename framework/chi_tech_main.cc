#include "chi_runtime.h"

//######################################################### Program entry point
/** Program entry point.

\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.

*/
int main(int argc, char** argv)
{
  chi::Initialize(argc, argv);

  int error_code;
  if (chi::run_time::sim_option_interactive_)
    error_code = chi::RunInteractive(argc, argv);
  else
    error_code = chi::RunBatch(argc, argv);

  chi::Finalize();

  return error_code;
}