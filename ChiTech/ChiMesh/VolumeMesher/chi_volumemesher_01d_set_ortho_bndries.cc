#include "ChiMesh/VolumeMesher/chi_volumemesher.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Sets boundary numbers on boundaries orthogonal to the cardinal directions
 * as xmax=0, xmin=1, ymax=2, ymin=3, zmax=4, zmin=5.*/
void chi_mesh::VolumeMesher::
SetupOrthogonalBoundaries()
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Setting orthogonal boundaries.";

  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  auto vol_cont = handler.GetGrid();

  const chi_mesh::Vector3 ihat(1.0, 0.0, 0.0);
  const chi_mesh::Vector3 jhat(0.0, 1.0, 0.0);
  const chi_mesh::Vector3 khat(0.0, 0.0, 1.0);

  for (auto& cell : vol_cont->local_cells)
    for (auto& face : cell.faces)
      if (not face.has_neighbor)
      {
        chi_mesh::Vector3& n = face.normal;

        int boundary_id = -1;
        if      (n.Dot(ihat)>0.999)  boundary_id = 0;
        else if (n.Dot(ihat)<-0.999) boundary_id = 1;
        else if (n.Dot(jhat)> 0.999) boundary_id = 2;
        else if (n.Dot(jhat)<-0.999) boundary_id = 3;
        else if (n.Dot(khat)> 0.999) boundary_id = 4;
        else if (n.Dot(khat)<-0.999) boundary_id = 5;

        if (boundary_id >= 0) face.neighbor_id = boundary_id;
      }//if bndry

  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Done setting orthogonal boundaries.";
}