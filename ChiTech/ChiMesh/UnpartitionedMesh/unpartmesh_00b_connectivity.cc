#include "chi_unpartitioned_mesh.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Establishes neighbor connectivity for the light-weight mesh.*/
void chi_mesh::UnpartitionedMesh::BuildMeshConnectivity()
{
  //======================================== Reset all cell neighbors
  int num_bndry_faces = 0;
  for (auto& cell : raw_cells)
    for (auto& face : cell->faces)
      if (not face.has_neighbor) ++num_bndry_faces;

  chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                              << " Number of unconnected faces "
                                 "before connectivity: " << num_bndry_faces;

  chi_log.Log() << "Establishing cell connectivity.";

  //======================================== Establish internal connectivity
  // Populate vertex subscriptions to internal cells
  vertex_cell_subscriptions.resize(vertices.size());
  uint64_t cur_cell_id=0;
  for (auto& cell : raw_cells)
  {
    for (auto vid : cell->vertex_ids)
      vertex_cell_subscriptions[vid].insert(cur_cell_id);
    ++cur_cell_id;
  }

  // Process raw cells
  std::set<size_t> cells_to_search; //This will be used and abused below
  cur_cell_id=0;
  for (auto& cell : raw_cells)
  {
    cells_to_search.clear();
    for (uint64_t vid : cell->vertex_ids)
      for (uint64_t cell_id : vertex_cell_subscriptions[vid])
        if (cell_id != cur_cell_id)
          cells_to_search.insert(cell_id);

    for (auto& cur_cell_face : cell->faces)
    {
      if (cur_cell_face.has_neighbor) continue;

      std::set<uint64_t> cfvids(cur_cell_face.vertex_ids.begin(),
                                cur_cell_face.vertex_ids.end());

      for (uint64_t adj_cell_id : cells_to_search)
      {
        auto adj_cell = raw_cells[adj_cell_id];

        for (auto& adj_cell_face : adj_cell->faces)
        {
          if (adj_cell_face.has_neighbor) continue;
          std::set<uint64_t> afvids(adj_cell_face.vertex_ids.begin(),
                                    adj_cell_face.vertex_ids.end());

          if (cfvids == afvids)
          {
            cur_cell_face.neighbor = adj_cell_id;
            adj_cell_face.neighbor = cur_cell_id;

            cur_cell_face.has_neighbor = true;
            adj_cell_face.has_neighbor = true;

            goto face_neighbor_found;
          }
        }//for adjacent cell face
      }
      face_neighbor_found:;
    }//for face

    ++cur_cell_id;
  }//for cell

  chi_log.Log() << "Establishing cell boundary connectivity.";

  //======================================== Establish boundary connectivity
  // Make list of internal cells on the boundary
  std::vector<LightWeightCell*> internal_cells_on_boundary;
  for (auto& cell : raw_cells)
  {
    bool cell_on_boundary = false;
    for (auto& face : cell->faces)
      if (not face.has_neighbor)
      { cell_on_boundary = true; break; }

    if (cell_on_boundary) internal_cells_on_boundary.push_back(cell);
  }

  // Populate vertex subscriptions to boundary cells
  std::vector<std::set<uint64_t>>
    vertex_bndry_cell_subscriptions(vertices.size());
  cur_cell_id=0;
  for (auto& cell : raw_boundary_cells)
  {
    for (auto vid : cell->vertex_ids)
      vertex_bndry_cell_subscriptions[vid].insert(cur_cell_id);
    ++cur_cell_id;
  }

  // Process boundary cells
  cur_cell_id=0;
  for (auto& cell : internal_cells_on_boundary)
    for (auto& face : cell->faces)
    {
      if (face.has_neighbor) continue;
      std::set<uint64_t> cfvids(face.vertex_ids.begin(),
                                face.vertex_ids.end());

      cells_to_search.clear();
      for (uint64_t vid : face.vertex_ids)
        for (uint64_t cell_id : vertex_bndry_cell_subscriptions[vid])
          cells_to_search.insert(cell_id);

      for (uint64_t adj_cell_id : cells_to_search)
      {
        auto& adj_cell = raw_boundary_cells[adj_cell_id];

        std::set<uint64_t> afvids(adj_cell->vertex_ids.begin(),
                                  adj_cell->vertex_ids.end());

        if (cfvids == afvids)
        {
          face.neighbor = adj_cell->material_id;
          break;
        }
      }//for adj_cell_id
    }//for face

  chi_log.Log() << "Done establishing cell connectivity.";

  num_bndry_faces = 0;
  for (auto cell : raw_cells)
    for (auto& face : cell->faces)
      if (not face.has_neighbor) ++num_bndry_faces;

  chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                              << " Number of boundary faces "
                                 "after connectivity: " << num_bndry_faces;
}