#include"ChiLua/chi_lua.h"

#include"../Solver/diffusion_solver.h"
#include "../Boundaries/chi_diffusion_bndry_reflecting.h"
#include "../Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "../Boundaries/chi_diffusion_bndry_robin.h"

#include"ChiPhysics/chi_physics.h"
#include <chi_log.h>

extern ChiPhysics&  chi_physics_handler;
extern ChiLog& chi_log;

#include "ChiMath/SpatialDiscretization/PiecewiseLinear/pwl.h"

#define DISCRETIZATION_METHOD 1
#define MAX_ITERS             2
#define RESIDUAL_TOL          3
#define BOUNDARY_TYPE         4
#define PROPERTY_D_MAP        5
#define PROPERTY_Q_MAP        6
#define PROPERTY_SIGMAA_MAP   7

//#############################################################################
/** Sets a property of a Diffusion solver. Please also consult the whitepaper
 * for the Diffusion solver (<a
 * href="../../whitepages/DiffusionSolver/DiffusionSolver.pdf">
 * Diffusion Whitepaper</a>)



\param SolverHandle int Handle to an existing diffusion solver.
\param PropertyIndex int Code for a specific property.
\param Values varying Number of inputs associated with the index.<br>

##_

###PropertyIndex\n
DISCRETIZATION_METHOD\n
 Discretization method. Expects to be followed by <B>DiscretizationMethod</B>
  (see below).\n\n

MAX_ITERS\n
 Solver maximum number of iterations.\n\n

RESIDUAL_TOL\n
 Residual tolerance. Expects to be followed by a floating point
 value (i.e. 1.0e-6).\n\n

BOUNDARY_TYPE\n
 Boundary type. Expects boundary index then <B>BoundaryTypeIndex</B>
 then type value.\n\n


PROPERTY_D_MAP\n
 Followed by an integer, sets the property index of a material
 from where the solver will collect the diffusion coefficient.\n\n

PROPERTY_Q_MAP\n
 Followed by an integer, sets the property index of a material
 from where the solver will collect the volumetric source value.\n\n

PROPERTY_SIGMAA_MAP\n
 Followed by an integer, sets the property index of a material
 from where the solver will collect the absorbtion cross-section
                 \f$ \sigma_a \f$.\n\n

\code
chiDiffusionSetProperty(solver,BOUNDARY_TYPE,2,DIFFUSION_DIRICHLET,1.0)
\endcode

###DiscretizationMethod\n
 PWLC\n
 Piecewise Linear Finite Element Continuous.\n\n

 PWLD_MIP\n
 Piecewise Linear Finite Element Discontinuous using the Modified
            Interior Penalty (MIP) method.\n\n


### BoundaryTypeIndex
DIFFUSION_REFLECTING
 Reflecting boundary conditions. Synonymous with Neumann with a
 derivative of 0.0.
             \f[ -D \hat{n}\cdot \nabla \phi = 0 \f]\n\n

DIFFUSION_DIRICHLET\n
 Constant value boundary condition.
 Expects to be followed by a value \f$ f \f$ associated with \f$ \phi \f$.
            \f[ \phi = f \f]\n\n

DIFFUSION_NEUMANN\n
 Constant derivative boundary condition. Expects to be followed
 by a constant \f$ f \f$ representing
                    \f[ -D \hat{n}\cdot \nabla \phi = f \f]\n\n

DIFFUSION_VACUUM\n
 Vacuum boundary conditions. More appropriate to neutron diffusion.
   \f[ \frac{1}{4}\phi + \frac{1}{2} D \hat{n}\cdot \nabla \phi = 0 \f]\n\n

DIFFUSION_ROBIN\n
 Robin boundary condition of the form
                   \f[ a \phi + b D \hat{n}\cdot \nabla \phi = f \f]\n\n

\ingroup LuaDiffusion
\author Jan*/
int chiDiffusionSetProperty(lua_State *L)
{
  int num_args = lua_gettop(L);

  //============================================= Get solver
  int solver_index = lua_tonumber(L,1);
  chi_diffusion::Solver* solver;

  try{
    solver = (chi_diffusion::Solver*)chi_physics_handler.solver_stack.at(solver_index);
  }
  catch(const std::out_of_range& o){
    chi_log.Log(LOG_0ERROR) << "ERROR: Invalid solver handle." << std::endl;
    exit(EXIT_FAILURE);
  }

  //============================================= Get property index
  int property = lua_tonumber(L,2);

  //============================================= Handle properties
  if (property == DISCRETIZATION_METHOD)
  {
    int method = lua_tonumber(L,3);
    if (method == PWLC)
    {
      SpatialDiscretization_PWL* discretization = new SpatialDiscretization_PWL;
      solver->discretization = discretization;
      solver->fem_method = PWLC;
    }
    else if (method == PWLD_MIP)
    {
      SpatialDiscretization_PWL* discretization = new SpatialDiscretization_PWL;
      solver->discretization = discretization;
      solver->fem_method = PWLD_MIP;
    }
    else
    {
      chi_log.Log(LOG_0ERROR) << "Invalid option for Discretization method in "
                   "chiDiffusionSetProperty.";
      exit(EXIT_FAILURE);
    }
  }
  else if (property == MAX_ITERS)
  {
    int num_iters = lua_tonumber(L,3);
    solver->max_iters = num_iters;
    chi_log.Log(LOG_ALLVERBOSE_2)
      << "Diffusion Solver: Maximum iterations set to " << solver->max_iters;
  }
  else if (property == RESIDUAL_TOL)
  {
    double tol = lua_tonumber(L,3);
    solver->residual_tolerance = tol;
  }

  else if (property == BOUNDARY_TYPE)
  {
    if (not solver->common_items_initialized)
      solver->InitializeCommonItems();

    chi_log.Log(LOG_0) << "Number of accessable boundaries: "
                       << solver->boundaries.size();

    if (num_args < 4)
    {
      chi_log.Log(LOG_0ERROR)
      << "Invalid amount of arguments used in"
      << " chiDiffusionSetProperty(...,BOUNDARY_TYPE.... "
      << " At least 4 arguments are expected.";
      exit(EXIT_FAILURE);
    }

    int bound_index = lua_tonumber(L,3);
    if (bound_index >= solver->boundaries.size())
    {
      chi_log.Log(LOG_0ERROR)
      << "Invalid boundary handle used in"
      << " chiDiffusionSetProperty(...,BOUNDARY_TYPE....";
      exit(EXIT_FAILURE);
    }

    int type_index = lua_tonumber(L,4);

    if (type_index == DIFFUSION_REFLECTING)
    {
      if (num_args != 4)
      {
        chi_log.Log(LOG_0ERROR)
          << "Invalid amount of arguments used in"
          << " chiDiffusionSetProperty(...,BOUNDARY_TYPE,"
          << bound_index << ",DIFFUSION_REFLECTING. "
          << " 4 arguments are expected.";
        exit(EXIT_FAILURE);
      }

      chi_diffusion::BoundaryReflecting* bound =
        new chi_diffusion::BoundaryReflecting;

      solver->boundaries[bound_index] = bound;

      chi_log.Log(LOG_0) << "Boundary " << bound_index << " set as "
                         << "Reflecting.";
    }
    else if (type_index == DIFFUSION_DIRICHLET)
    {
      if (num_args != 5)
      {
        chi_log.Log(LOG_0ERROR)
          << "Invalid amount of arguments used in"
          << " chiDiffusionSetProperty(...,BOUNDARY_TYPE,"
          << bound_index << ",DIFFUSION_DIRICHLET. "
          << " 5 arguments are expected.";
        exit(EXIT_FAILURE);
      }

      double b_value = lua_tonumber(L,5);

      chi_diffusion::BoundaryDirichlet* bound =
        new chi_diffusion::BoundaryDirichlet;
      bound->boundary_value = b_value;

      solver->boundaries[bound_index] = bound;

      chi_log.Log(LOG_0) << "Boundary " << bound_index << " set as "
                         << "Dirichlet with value " << b_value;
    }
    else if (type_index == DIFFUSION_NEUMANN)
    {
      if (num_args != 5)
      {
        chi_log.Log(LOG_0ERROR)
          << "Invalid amount of arguments used in"
          << " chiDiffusionSetProperty(...,BOUNDARY_TYPE,"
          << bound_index << ",DIFFUSION_NEUMANN. "
          << " 5 arguments are expected.";
        exit(EXIT_FAILURE);
      }

      double f_value = lua_tonumber(L,5);

      chi_diffusion::BoundaryRobin* bound =
        new chi_diffusion::BoundaryRobin(0.0,1.0,f_value);

      solver->boundaries[bound_index] = bound;

      chi_log.Log(LOG_0) << "Boundary " << bound_index << " set as "
                         << "Neumann with f = ("
                         << f_value << ") ";
    }
    else if (type_index == DIFFUSION_VACUUM)
    {
      if (num_args != 4)
      {
        chi_log.Log(LOG_0ERROR)
          << "Invalid amount of arguments used in"
          << " chiDiffusionSetProperty(...,BOUNDARY_TYPE,"
          << bound_index << ",DIFFUSION_VACUUM. "
          << " 4 arguments are expected.";
        exit(EXIT_FAILURE);
      }

      chi_diffusion::BoundaryRobin* bound =
        new chi_diffusion::BoundaryRobin(0.25,0.5,0.0);

      solver->boundaries[bound_index] = bound;

      chi_log.Log(LOG_0) << "Boundary " << bound_index << " set as "
                         << "Vacuum.";
    }
    else if (type_index == DIFFUSION_ROBIN)
    {
      if (num_args != 7)
      {
        chi_log.Log(LOG_0ERROR)
          << "Invalid amount of arguments used in"
          << " chiDiffusionSetProperty(...,BOUNDARY_TYPE,"
          << bound_index << ",DIFFUSION_ROBIN. "
          << " 7 arguments are expected.";
        exit(EXIT_FAILURE);
      }

      double a_value = lua_tonumber(L,5);
      double b_value = lua_tonumber(L,6);
      double f_value = lua_tonumber(L,7);

      chi_diffusion::BoundaryRobin* bound =
        new chi_diffusion::BoundaryRobin(a_value,b_value,f_value);

      solver->boundaries[bound_index] = bound;

      chi_log.Log(LOG_0) << "Boundary " << bound_index << " set as "
                         << "Robin with a,b,f = ("
                         << a_value << ","
                         << b_value << ","
                         << f_value << ") ";
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Unsupported boundary type encountered in call to "
        << "chiDiffusionSetProperty(..,BOUNDARY_TYPE,.. :"
        << type_index;
      exit(EXIT_FAILURE);
    }
  }
  else if (property == PROPERTY_D_MAP)
  {
    solver->property_map_D = lua_tonumber(L,3);
  }
  else if (property == PROPERTY_Q_MAP)
  {
    solver->property_map_q = lua_tonumber(L,3);
  }
  else if (property == PROPERTY_SIGMAA_MAP)
  {
    solver->property_map_sigma = lua_tonumber(L,3);
  }
  else
  {
    chi_log.Log(LOG_0ERROR) << "Invalid property in chiDiffusionSetProperty.";
    exit(EXIT_FAILURE);
  }
  return 0;
}
