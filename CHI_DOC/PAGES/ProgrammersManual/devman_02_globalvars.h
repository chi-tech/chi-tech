/**\page DevManGlobalVars Global variables

\section devman2_sec0 Global variables available

The major global variables are defined in "chi_tech_main.h"

They are:
 - CHI_CONSOLE chi_console The link to the lua scripting engine.
 - CHI_MATH chi_math_handler A handler of math related entities
 - CHI_PHYSICS chi_physics_handler A handler of physics related items
 - CHI_MPI chi_mpi A handler for parallel related items.
 - CHI_LOG chi_log A handler for parallel logging events and verbosity.




\subsection devman2_sec0_3 Connecting to the physics handler

The physics handler maintains a number of data structures, most notably
 are the following three:
 - CHI_PHYSICS::solver_stack Solvers are pushed here
 - CHI_PHYSICS::material_stack Materials are pushed here
 - CHI_PHYSICS::fieldfunc_stack Field functions are pushed here

To access chi_physics_handler include the following code at the top of your
 code

\code
#include <CHI_PHYSICS/chi_physics.h>

extern CHI_PHYSICS chi_physics_handler;
\endcode



\subsection devman2_sec0_4 Connecting to MPI

General MPI information like the current location id and the total amount
 of parallel processes is contained in CHI_MPI:
 - CHI_MPI::location_id
 - CHI_MPI::process_count

Additionally, by including the headers for chi_mpi, developers have access to
 all the classic mpi headers.

\code
#include <chi_mpi.h>

extern CHI_MPI chi_mpi;

if (chi_mpi.location_id == 1)
 printf("This is process 1, Dab!");
\endcode



\subsection devman2_sec0_5 Connecting to the parallel logging utility

Printing information in a parallel environment can be a very involved
process. One can't simply use "std::cout <<" on every process otherwise
the output to the log will be chaotic. For this reason we employ a common
 logging utility which return an output string stream using the function
 call CHI_LOG::Log.

Connecting to chi_log is done as follows
\code
#include <chi_log.h>

extern CHI_LOG chi_log;
\endcode

The logger needs to be supplied with an enumeration (LOG_LVL) indicating
 the type of output. The following enumerations are supported:
 - LOG_0,                      Used only for location 0
 - LOG_0WARNING,               Warning only for location 0
 - LOG_0ERROR,                 Error only for location 0
 - LOG_0VERBOSE_0,             Default verbosity level
 - LOG_0VERBOSE_1,             Used only if verbosity level equals 1
 - LOG_0VERBOSE_2,             Used only if verbosity level equals 2
 - LOG_ALL,                    Verbose level 0 all locations
 - LOG_ALLWARNING,             Warning for any location
 - LOG_ALLERROR,               Error for any location
 - LOG_ALLVERBOSE_0,     Default verbosity level
 - LOG_ALLVERBOSE_1,     Used only if verbosity level equals 1
 - LOG_ALLVERBOSE_2,     Used only if verbosity level equals 2



\todo Commandline parameter documentation


*/