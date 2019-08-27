#include "lbs_linear_boltzman_solver.h"
#include <CHI_MESH/CHI_MESHHANDLER/chi_meshhandler.h>
#include <CHI_MESH/CHI_VOLUMEMESHER/chi_volumemesher.h>
#include <CHI_MESH/CHI_VOLUMEMESHER/Linemesh1D/volmesher_linemesh1d.h>
#include <CHI_MESH/CHI_VOLUMEMESHER/Extruder/volmesher_extruder.h>
#include <CHI_MESH/CHI_VOLUMEMESHER/Predefined2D/volmesher_predefined2d.h>


#include <chi_mpi.h>
#include <chi_log.h>

extern CHI_MPI chi_mpi;
extern CHI_LOG chi_log;

typedef chi_mesh::SweepManagement::AngleSet TAngleSet;
typedef chi_mesh::SweepManagement::AngleSetGroup TAngleSetGroup;

#include <iomanip>
#include "../../CHI_CONSOLE/chi_console.h"

extern CHI_CONSOLE chi_console;

//###################################################################
/**Initializes the sweep ordering for the given groupset.*/
void CHI_NPTRANSPORT::ComputeSweepOrderings(NPT_GROUPSET *groupset)
{
  chi_log.Log(LOG_0) << "Computing Sweep ordering.\n";

  //============================================= Clear sweep ordering
  sweep_orderings.clear();
  sweep_orderings.shrink_to_fit();

  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D MESHES
  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
  {
    int num_azi = groupset->quadrature->azimu_ang.size();
    int num_pol = groupset->quadrature->polar_ang.size();
    int pa      = num_pol/2;

    if (num_azi != 1)
    {
      chi_log.Log(LOG_0)
        << "Incompatible number of azimuthal angles in quadrature set "
        << "for a 1D simulation.";
      exit(EXIT_FAILURE);
    }

    chi_mesh::SweepManagement::SPDS* new_swp_order =
      chi_mesh::SweepManagement::
      CreateSweepOrder(groupset->quadrature->polar_ang[0],
                       groupset->quadrature->azimu_ang[0],
                       this->grid,
                       groupset->groups.size());
    this->sweep_orderings.push_back(new_swp_order);

    new_swp_order =
      chi_mesh::SweepManagement::
      CreateSweepOrder(groupset->quadrature->polar_ang[pa],
                       groupset->quadrature->azimu_ang[0],
                       this->grid,
                       groupset->groups.size());
    this->sweep_orderings.push_back(new_swp_order);
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D MESHES
  else if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherPredefined2D))
  {
    int num_azi = groupset->quadrature->azimu_ang.size();
    int num_pol = groupset->quadrature->polar_ang.size();

    if (num_pol != 1)
    {
      chi_log.Log(LOG_0)
        << "Incompatible number of polar angles in quadrature set "
        << "for a 2D simulation.";
      exit(EXIT_FAILURE);
    }

    for (int i=0; i<num_azi; i++)
    {
      chi_mesh::SweepManagement::SPDS* new_swp_order =
        chi_mesh::SweepManagement::
        CreateSweepOrder(groupset->quadrature->polar_ang[0],
                         groupset->quadrature->azimu_ang[i],
                         this->grid,
                         groupset->groups.size());
      this->sweep_orderings.push_back(new_swp_order);
    }
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRUDED MESHES
  else if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherExtruder))
  {
    int num_azi = groupset->quadrature->azimu_ang.size();
    int num_pol = groupset->quadrature->polar_ang.size();

    if (num_azi < 4)
    {
      chi_log.Log(LOG_0)
        << "Incompatible number of azimuthal angles in quadrature set "
        << "for a 3D simulation.";
      exit(EXIT_FAILURE);
    }
    if (num_pol < 2)
    {
      chi_log.Log(LOG_0)
        << "Incompatible number of polar angles in quadrature set "
        << "for a 3D simulation.";
      exit(EXIT_FAILURE);
    }

    //============================================= Create sweep ordering
    //                                              per azimuthal angle
    //                                              per hemisphere

    //=========================================== TOP HEMISPHERE
    for (int i=0; i<num_azi; i++)
    {
      chi_mesh::SweepManagement::SPDS* new_swp_order =
        chi_mesh::SweepManagement::
                  CreateSweepOrder(groupset->quadrature->polar_ang[0],
                                   groupset->quadrature->azimu_ang[i],
                                   this->grid,
                                   groupset->groups.size());
      this->sweep_orderings.push_back(new_swp_order);
    }
    //=========================================== BOTTOM HEMISPHERE
    int pa = num_pol/2;
    for (int i=0; i<num_azi; i++)
    {
      chi_mesh::SweepManagement::SPDS* new_swp_order =
        chi_mesh::SweepManagement::
        CreateSweepOrder(groupset->quadrature->polar_ang[pa],
                         groupset->quadrature->azimu_ang[i],
                         this->grid,
                         groupset->groups.size());
      this->sweep_orderings.push_back(new_swp_order);
    }

  }
  else
  {
    fprintf(stderr,"ERROR: Cannot create sweep ordering"
                   " for given mesh type.\n");
    exit(EXIT_FAILURE);
  }


  chi_log.Log(LOG_0)
    << "Done computing sweep orderings.           Process memory = "
    << std::setprecision(3)
    << chi_console.GetMemoryUsageInMB() << " MB";

}
