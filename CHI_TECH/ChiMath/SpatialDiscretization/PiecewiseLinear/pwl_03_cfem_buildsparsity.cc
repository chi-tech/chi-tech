#include "pwl.h"

#include "chi_log.h"
extern ChiLog chi_log;

#include "chi_mpi.h"
extern ChiMPI chi_mpi;

//###################################################################
/**Builds the sparsity pattern for a Continuous Finite Element Method.*/
void SpatialDiscretization_PWL::
  BuildCFEMSparsityPattern(chi_mesh::MeshContinuum *grid,
  std::vector<int> &nodal_bndry_ids,
  std::vector<int> &nodal_nnz_in_diag,
  std::vector<int> &nodal_nnz_off_diag,
  const std::pair<int,int>& domain_ownership)
{
  std::vector<std::vector<int>> nodal_connections(grid->nodes.size());

  nodal_bndry_ids   .resize(grid->nodes.size(),0);
  nodal_nnz_in_diag .resize(grid->nodes.size(),0);
  nodal_nnz_off_diag.resize(grid->nodes.size(),0);

  for (auto& glob_index : grid->local_cell_glob_indices)
  {
    auto cell = grid->cells[glob_index];

    for (size_t i=0; i < cell->vertex_ids.size(); i++)
    {
      int ir =  MapNode(cell->vertex_ids[i]);

      if (ir<0)
      {
        chi_log.Log(LOG_ALLERROR)
          << "ir Mapping error node " << cell->vertex_ids[i];
        exit(EXIT_FAILURE);
      }

      //================================== Check if i is on boundary
      for (auto& face : cell->faces)
      {
        if (face.neighbor < 0)
        {
          for (auto& fvid : face.vertex_ids)
          {
            int v0_index = MapNode(fvid);

            if (v0_index<0)
            {
              chi_log.Log(LOG_ALLERROR)
                << "v0 Mapping error node " << fvid;
              exit(EXIT_FAILURE);
            }

            if (ir == v0_index)
            {
              nodal_bndry_ids[ir]= face.neighbor;
              break;
            } //if ir part of face
          }
        }
      }//for f

      //======================================= Set nodal connections
      std::vector<int>& node_links = nodal_connections[ir];
      for (auto& jvid : cell->vertex_ids)
      {
        int jr = MapNode(jvid);

        //====================== Check for duplicates
        bool already_there = false;
        for (auto& node_link : node_links)
        {
          if (node_link == jr)
          {already_there = true; break;}
        }
        if (!already_there)
        {
          node_links.push_back(jr);
          if ((jr>=domain_ownership.first) and
              (jr<=domain_ownership.second))
          {
            nodal_nnz_in_diag[ir]+=1;
          } else
          {
            nodal_nnz_off_diag[ir]+=1;
          }
        }
      }//for j
    }//for i
  }

  MPI_Barrier(MPI_COMM_WORLD);
}