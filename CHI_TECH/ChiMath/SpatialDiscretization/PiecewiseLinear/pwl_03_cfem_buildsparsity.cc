#include "pwl.h"

#include "chi_log.h"
extern ChiLog chi_log;

#include "chi_mpi.h"
extern ChiMPI chi_mpi;

//###################################################################
/**Builds the sparsity pattern for a Continuous Finite Element Method.*/
void SpatialDiscretization_PWL::
  BuildCFEMSparsityPattern(chi_mesh::MeshContinuum *grid,
                           std::vector<int> &nodal_nnz_in_diag,
                           std::vector<int> &nodal_nnz_off_diag,
                           const std::pair<int,int>& domain_ownership)
{
  std::vector<std::vector<int>> nodal_connections(grid->vertices.size());

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag .resize(grid->vertices.size(), 0);
  nodal_nnz_off_diag.resize(grid->vertices.size(), 0);

  for (auto& cell : grid->local_cells)
  {
    for (size_t i=0; i < cell.vertex_ids.size(); i++)
    {
      int ir = MapCFEMDOF(cell.vertex_ids[i]);

      if (ir<0)
      {
        chi_log.Log(LOG_ALLERROR)
          << "ir Mapping error node " << cell.vertex_ids[i];
        exit(EXIT_FAILURE);
      }

      //======================================= Set nodal connections
      std::vector<int>& node_links = nodal_connections[ir];
      for (auto& jvid : cell.vertex_ids)
      {
        int jr = MapCFEMDOF(jvid);

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
          if ((jr>=cfem_local_block_address) and
              (jr<=(cfem_local_block_address+domain_ownership.first)))
          {
            nodal_nnz_in_diag[ir]+=1;
          } else
          {
            nodal_nnz_off_diag[ir]+=1;
          }
        }
      }//for j
    }//for i
  }//for cell

  MPI_Barrier(MPI_COMM_WORLD);
}