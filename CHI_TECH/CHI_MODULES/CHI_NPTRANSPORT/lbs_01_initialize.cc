#include "linear_boltzman_solver.h"
#include "../../CHI_MESH/CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../CHI_MESH/CHI_CELL/cell.h"
#include "../../CHI_MATH/CHI_DISCRETIZATION/CHI_DISCRETIZATION_PWL/pwl.h"


#include <chi_mpi.h>
#include <chi_log.h>

extern CHI_MPI chi_mpi;
extern CHI_LOG chi_log;

#include <iomanip>
#include "../../CHI_CONSOLE/chi_console.h"

extern CHI_CONSOLE chi_console;



//###################################################################
/** Initialize the solver.*/
void CHI_NPTRANSPORT::Initialize()
{
  //============================================= Input checks
  if (discretization == nullptr)
  {
    chi_log.Log(LOG_ALLERROR)
      << "CHI_NPTRANSPORT: No discretization method set.";
    exit(EXIT_FAILURE);
  }

//  int L = options.scattering_order;
//  this->num_moments = L*(L+2) + 1;
  ComputeNumberOfMoments();

  if (chi_mpi.location_id == 0)
  {
    chi_log.Log(LOG_0) << "Initializing NPT Solver\n";
    chi_log.Log(LOG_0) << "Scattering order    : "
                       << options.scattering_order << std::endl;
    chi_log.Log(LOG_0) << "Number of Groups    : "
                       << groups.size() << std::endl;
    chi_log.Log(LOG_0) << "Number of Group sets: "
                       << group_sets.size() << std::endl;
    chi_log.Log(LOG_0) << "Number of Moments   : "
                       << num_moments << std::endl;



    //================================================== Output Groupsets
    for (int gs=0; gs<group_sets.size(); gs++)
    {
      char buf_pol[10];
      std::string outstr;
      int counter = 0;

      chi_log.Log(LOG_0) << "\n***** Groupset " << gs << " *****\n";
      chi_log.Log(LOG_0) << "Groups: ";
      outstr = std::string("");
      counter = 0;
      for (int g=0; g<group_sets[gs]->groups.size(); g++)
      {
        snprintf(buf_pol,10,"%5d ",group_sets[gs]->groups[g]->id);
        outstr += std::string(buf_pol);
        counter++;
        if (counter == 12)
        {
          counter = 0;
          chi_log.Log(LOG_0) << outstr << "\n";
          outstr = std::string("");
        }

      }//for g
      chi_log.Log(LOG_0) << outstr << "\n\n";


      CHI_PRODUCT_QUADRATURE* quad = group_sets[gs]->quadrature;
      chi_log.Log(LOG_0VERBOSE_1) << "Product Quadrature polar angles:\n";
      outstr = std::string("");
      counter = 0;
      for (int ang=0; ang<quad->polar_ang.size(); ang++)
      {
        snprintf(buf_pol,10,"%9.3f ",quad->polar_ang[ang]);
        outstr += std::string(buf_pol);
        counter++;
        if (counter == 6)
        {
          counter = 0;
          chi_log.Log(LOG_0VERBOSE_1) << outstr << "\n";
          outstr = std::string("");
        }
      }
      chi_log.Log(LOG_0VERBOSE_1) << outstr << "\n\n";


      chi_log.Log(LOG_0VERBOSE_1) << "Product Quadrature azimuthal angles:\n";
      outstr = std::string("");
      counter = 0;
      for (int ang=0; ang<quad->azimu_ang.size(); ang++)
      {
        snprintf(buf_pol,10,"%9.3f ",quad->azimu_ang[ang]);
        outstr += std::string(buf_pol);
        counter++;
        if (counter == 6)
        {
          counter = 0;
          chi_log.Log(LOG_0VERBOSE_1) << outstr << "\n";
          outstr = std::string("");
        }
      }
      chi_log.Log(LOG_0VERBOSE_1) << outstr << "\n\n";
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //================================================== Determine partitioning
  //SetPartitioning();


  //================================================== Compute cell matrices
  chi_log.Log(LOG_0) << "Computing cell matrices.\n";
  chi_mesh::Region*  aregion = this->regions.back();
  this->grid                 = aregion->volume_mesh_continua.back();


  discretization->AddViewOfLocalContinuum(grid,
                                          grid->local_cell_glob_indices.size(),
                                          grid->local_cell_glob_indices.data());

  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << "Cell matrices computed.                   Process memory = "
    << std::setprecision(3)
    << chi_console.GetMemoryUsageInMB() << " MB";


  //================================================== Add transport views
  CHI_DISCRETIZATION_PWL* pwl_discretization =
    (CHI_DISCRETIZATION_PWL*)discretization;

  std::set<int> unique_material_ids;
  int invalid_mat_cell_count = 0;
  for (int c=0; c<grid->local_cell_glob_indices.size(); c++)
  {
    int cell_g_index = grid->local_cell_glob_indices[c];
    auto cell = grid->cells[cell_g_index];

    //Transport simulations can sometimes have very complex material
    //specifications. Sometimes users pre-emptively add many more
    //materials than are required. We therefore here initialize
    //only materials used in the simulation.
    unique_material_ids.insert(cell->material_id);
    if (cell->material_id<0)
      invalid_mat_cell_count++;

    //Transport views act as a data structure to store information
    //related to the transport simulation. The most prominent function
    //here is that it holds the means to know where a given cell's
    //transport quantities are located in the unknown vectors (i.e. phi)
    CellFEView* cell_fe_view = pwl_discretization->MapFeView(cell_g_index);
    NPT_CELLVIEW_FULL* full_cell_view =
      new NPT_CELLVIEW_FULL(cell_fe_view->dofs, groups.size(), num_moments);
    cell_transport_views.push_back(full_cell_view);

    //For PWLD, for a given cell, within a given sweep chunk,
    // we need to solve a matrix which square size is the
    // amount of dofs on the cell. max_cell_dof_count is
    // initialized here.
    if (cell_fe_view->dofs > max_cell_dof_count)
      max_cell_dof_count = cell_fe_view->dofs;

  }

  if (invalid_mat_cell_count>0)
  {
    chi_log.Log(LOG_ALLWARNING)
      << "Number of invalid material cells: " << invalid_mat_cell_count;
  }



  //================================================== Initialize materials
  InitMaterials(unique_material_ids);


  //================================================== Initialize parrays
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << "Initializing parallel arrays. " << std::endl;

  InitializeParrays();


  chi_log.Log(LOG_0)
    << "Done with parallel arrays.                Process memory = "
    << std::setprecision(3)
    << chi_console.GetMemoryUsageInMB() << " MB" << std::endl;

  //================================================== Initialize communicators
  InitializeCommunicators();

  MPI_Barrier(MPI_COMM_WORLD);
}