#ifndef _chi_tech_main_h
#define _chi_tech_main_h


#include"CHI_CONSOLE/chi_console.h"
CHI_CONSOLE       chi_console;
CHI_VECTOR<char>  chiconsoleInputBuffer;

#ifdef CHI_USEGRAPHICS
  #include <omp.h>
  #include "CHI_GRAPHICS/chi_graphics.h"
  #include "CHI_WINDOWMANAGER/chi_windowmanager.h"
  CHI_GRAPHICS chigraphics;
  CHI_WINDOWMANAGER chiwindowManager;

  #include"CHI_GRAPHICS/CHI_MATERIAL/chi_material.h"
  #include"CHI_GRAPHICS/CHI_OBJECT/chi_object.h"
  #include"CHI_GRAPHICS/CHI_SURFACE/chi_surface.h"
  CHI_VECTOR<CHI_MATERIAL>        chimaterialStack;
  CHI_VECTOR<CHI_OBJECT>          chiobjectStack;
  CHI_VECTOR<CHI_SURFACE>         chisurfaceMeshStack;
#endif

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





//CHI_VECTOR<CHI_TRANSFORM>       chitoolstransformStack;
#ifdef CHI_USESURFACEREMESHER
    #ifdef CHI_USEGRAPHICS
      #include"CHI_MODULES/CHI_SURFACEREMESHER/chi_surfaceremesher.h"
      CHI_VECTOR<CHI_SURFACEREMESHER> chisurfaceRemesherStack;
    #endif

#endif

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