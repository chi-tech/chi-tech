#include "volmesher_predefunpart.h"

#include "chi_log.h"

extern ChiLog& chi_log;

//###################################################################
/**Establishes neighbor connectivity for the light-weight mesh.*/
void chi_mesh::VolumeMesherPredefinedUnpartitioned::
  BuildMeshConnectivity(chi_mesh::UnpartitionedMesh* umesh)
{
  //======================================== Populate vertex
  //                                                   subscriptionns
  std::vector<std::set<int>> vertex_subs(umesh->vertices.size());
  int c=-1;
  for (auto cell : umesh->raw_cells)
  {
    ++c;
    for (auto vid : cell->vertex_ids)
      vertex_subs[vid].insert(c);
  }

  int num_bndry_faces = 0;
  for (auto cell : umesh->raw_cells)
    for (auto& face : cell->faces)
      if (face.neighbor < 0) ++num_bndry_faces;

  chi_log.Log(LOG_0) << "Number of bndry faces: " << num_bndry_faces;

  //======================================== Establish connectivity
  std::vector<std::set<int>> cells_to_search;
  c=-1;
  for (auto cell : umesh->raw_cells)
  {
    ++c;
    for (auto& face : cell->faces)
    {
      if (face.neighbor >= 0) continue;

      cells_to_search.clear();
      for (auto cfvid : face.vertex_ids)
        cells_to_search.push_back(vertex_subs[cfvid]);

      std::set<int> cfvids(face.vertex_ids.begin(),
                           face.vertex_ids.end());

      bool stop_searching = false;
      for (auto& adj_cell_id_set : cells_to_search)
      {
        for (auto& adj_cell_id : adj_cell_id_set)
        {
          if (adj_cell_id == c) continue;
          auto adj_cell = umesh->raw_cells[adj_cell_id];

          bool adj_cell_matches = false;
          for (auto& aface : adj_cell->faces)
          {
            std::set<int> afvids(aface.vertex_ids.begin(),
                                 aface.vertex_ids.end());

            if (afvids == cfvids)
            {
              adj_cell_matches = true;
              break;
            }
          }

          if (adj_cell_matches)
          {
            face.neighbor = adj_cell_id;
            stop_searching = true;
          }

          if (stop_searching) break;
        }//cell id
        if (stop_searching) break;
      }//cell id set
    }//for face
  }//for cell

  num_bndry_faces = 0;
  for (auto cell : umesh->raw_cells)
    for (auto& face : cell->faces)
      if (face.neighbor < 0) ++num_bndry_faces;

  chi_log.Log(LOG_0) << "Number of bndry faces: " << num_bndry_faces;
}