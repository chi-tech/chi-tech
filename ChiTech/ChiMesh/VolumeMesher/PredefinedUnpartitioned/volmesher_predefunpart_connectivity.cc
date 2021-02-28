#include "volmesher_predefunpart.h"

#include "chi_log.h"

extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Establishes neighbor connectivity for the light-weight mesh.*/
void chi_mesh::VolumeMesherPredefinedUnpartitioned::
  BuildMeshConnectivity(chi_mesh::UnpartitionedMesh* umesh)
{
  //======================================== Populate vertex
  //                                                   subscriptionns
  std::vector<std::set<size_t>> vertex_subs(umesh->vertices.size());
  uint64_t cur_cell_id=0;
  for (auto& cell : umesh->raw_cells)
  {
    for (auto vid : cell->vertex_ids)
      vertex_subs[vid].insert(cur_cell_id);
    ++cur_cell_id;
  }

  int num_bndry_faces = 0;
  for (auto cell : umesh->raw_cells)
    for (auto& face : cell->faces)
    {
      if (face.neighbor < 0) ++num_bndry_faces;
      face.neighbor = -1;
    }

  chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                              << " Number of boundary faces "
                                 "before connectivity: " << num_bndry_faces;

  //======================================== Establish connectivity
  chi_log.Log() << "Establishing cell connectivity.";
  std::set<size_t> cells_to_search;
  cur_cell_id=0;
  for (auto& cell : umesh->raw_cells)
  {
    cells_to_search.clear();
    for (uint64_t vid : cell->vertex_ids)
      for (uint64_t cell_id : vertex_subs[vid])
        if (cell_id != cur_cell_id)
          cells_to_search.insert(cell_id);

    for (auto& cur_cell_face : cell->faces)
    {
      if (cur_cell_face.neighbor >= 0 ) continue;

      std::set<uint64_t> cfvids(cur_cell_face.vertex_ids.begin(),
                                cur_cell_face.vertex_ids.end());

      for (uint64_t adj_cell_id : cells_to_search)
      {
        auto adj_cell = umesh->raw_cells[adj_cell_id];

        for (auto& adj_cell_face : adj_cell->faces)
        {
          if (adj_cell_face.neighbor >= 0) continue;
          std::set<uint64_t> afvids(adj_cell_face.vertex_ids.begin(),
                                    adj_cell_face.vertex_ids.end());

          if (cfvids == afvids)
          {
            cur_cell_face.neighbor = adj_cell_id;
            adj_cell_face.neighbor = cur_cell_id;
            goto face_neighbor_found;
          }
        }//for adjacent cell face
      }
      face_neighbor_found:;
    }//for face

    ++cur_cell_id;
  }//for cell

  chi_log.Log() << "Done establishing cell connectivity.";

  num_bndry_faces = 0;
  for (auto cell : umesh->raw_cells)
    for (auto& face : cell->faces)
      if (face.neighbor < 0) ++num_bndry_faces;

  chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                              << " Number of boundary faces "
                                 "after connectivity: " << num_bndry_faces;
}