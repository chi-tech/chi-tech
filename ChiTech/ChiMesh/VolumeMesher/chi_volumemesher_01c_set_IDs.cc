#include "chi_volumemesher.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"




//###################################################################
/**Sets material id's using a logical volume.*/
void chi_mesh::VolumeMesher::
SetMatIDFromLogical(const chi_mesh::LogicalVolume& log_vol,bool sense, int mat_id)
{
  chi_log.Log(LOG_0)
    << chi::program_timer.GetTimeString()
    << " Setting material id from logical volume.";
  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  chi_mesh::MeshContinuumPtr vol_cont = handler.GetGrid();

  int num_cells_modified = 0;
  for (auto& cell : vol_cont->local_cells)
  {
    if (log_vol.Inside(cell.centroid) && sense){
      cell.material_id = mat_id;
      ++num_cells_modified;
    }
  }

  const auto& ghost_ids = vol_cont->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = vol_cont->cells[ghost_id];
    if (log_vol.Inside(cell.centroid) && sense)
      cell.material_id = mat_id;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << chi::program_timer.GetTimeString()
    << " Done setting material id from logical volume. "
    << "Number of cells modified = " << num_cells_modified << ".";
}

//###################################################################
/**Sets material id's using a logical volume.*/
void chi_mesh::VolumeMesher::
SetBndryIDFromLogical(const chi_mesh::LogicalVolume& log_vol,bool sense, int bndry_id)
{
  chi_log.Log(LOG_0)
    << chi::program_timer.GetTimeString()
    << " Setting boundary id from logical volume.";
  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  chi_mesh::MeshContinuumPtr vol_cont = handler.GetGrid();

  int num_faces_modified = 0;
  for (auto& cell : vol_cont->local_cells)
  {
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor) continue;
      if (log_vol.Inside(face.centroid) && sense){
        face.neighbor_id = abs(bndry_id);
        ++num_faces_modified;
      }
    }
  }

  chi_log.Log(LOG_0)
    << chi::program_timer.GetTimeString()
    << " Done setting boundary id from logical volume. "
    << "Number of faces modified = " << num_faces_modified << ".";
}

//###################################################################
/**Sets material id's for all cells to the specified material id.*/
void chi_mesh::VolumeMesher::
SetMatIDToAll(int mat_id)
{
  chi_log.Log(LOG_0)
    << chi::program_timer.GetTimeString()
    << " Setting material id " << mat_id << "to all cells.";

  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  auto vol_cont = handler.GetGrid();

  for (auto& cell : vol_cont->local_cells)
    cell.material_id = mat_id;

  const auto& ghost_ids = vol_cont->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
    vol_cont->cells[ghost_id].material_id = mat_id;

  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << chi::program_timer.GetTimeString()
    << " Done setting material id " << mat_id << " to all cells";
}