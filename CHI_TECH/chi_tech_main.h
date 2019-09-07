#ifndef _chi_tech_main_h
#define _chi_tech_main_h


#include"CHI_CONSOLE/chi_console.h"
CHI_CONSOLE       chi_console;

#include "CHI_MATH/chi_math.h"
#include "CHI_PHYSICS/chi_physics.h"
#include "CHI_TIMER/chi_timer.h"
CHI_MATH    chi_math_handler;
CHI_PHYSICS chi_physics_handler;
CHI_TIMER   chi_program_timer;

#include <chi_mpi.h>
#include <chi_log.h>
CHI_MPI chi_mpi;
CHI_LOG chi_log;


//=============================================== Stacks

#include"CHI_TOOLS/CHI_TRANSFORM/chi_transform.h"
#include"CHI_MESH/CHI_MESHHANDLER/chi_meshhandler.h"


std::vector<chi_mesh::MeshHandler*>  chi_meshhandler_stack;
int                                  chi_current_mesh_handler=-1;


//###################################################################

//=============================================== Global variables
bool            chi_termination_posted = false;
std::string     input_file_name;
bool            sim_option_withgraphics = true;

/** GLOBAL TIMING
   This array serves as a place holder for many things
   [ 0]      Unused
       ...
   [09]      Maximum memory
   [10]      Memory events
   [11]      Memory event counter
   [12],[13] [12]*8/[13] gives the avg sweep buffer
             communication requirement in bytes
   [14],[15] [14]/[15] average sweep predecessor
             message count per angleset
   [16]      Cumulative sweep time
   [17]      Number of sweeps counter
   [18]      Cumulative set-source time
   [19]      Number of set source counts
*/
double          chi_global_timings[20];


#endif