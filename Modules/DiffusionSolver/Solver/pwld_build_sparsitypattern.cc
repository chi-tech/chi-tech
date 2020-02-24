#include "diffusion_solver.h"

#include <chi_log.h>
extern ChiLog chi_log;

#include <chi_mpi.h>
extern ChiMPI chi_mpi;

//###################################################################
/** Builds the sparsity pattern for PWLC spatial discretization.*/
void chi_diffusion::Solver::PWLDBuildSparsityPattern()
{
  int num_loc_cells = grid->local_cell_glob_indices.size();
  int dof_count = 0;
  std::set<int> local_border_cells;
  for (auto& cell : grid->local_cells)
  {
    auto ip_view = new DiffusionIPCellView;
    ip_view->cell_dof_start = dof_count + pwld_local_dof_start;
    pwld_cell_dof_array_address.push_back(dof_count);
    ip_cell_views.push_back(ip_view);

    for (size_t v=0; v<cell.vertex_ids.size(); v++)
    {
      nodal_nnz_in_diag[dof_count] = cell.vertex_ids.size();

      for (size_t f=0; f<cell.faces.size(); f++)
      {
        if (cell.faces[f].neighbor >= 0) //Not bndry
        {
          bool is_local = cell.faces[f].IsNeighborLocal(grid);

          if (is_local)
          {
            int neighbor_local_id = cell.faces[f].GetNeighborLocalID(grid);
            auto adj_cell = grid->local_cells[neighbor_local_id];
            nodal_nnz_in_diag[dof_count] += adj_cell.vertex_ids.size();
          }
          else
          {
            local_border_cells.insert(cell.local_id);
          }
        }
      }
      dof_count++;
    }

    //==================================== Boundary numbers
    for (size_t f=0; f<cell.faces.size(); f++)
    {
      if (cell.faces[f].neighbor < 0)
      {
        for (size_t fv=0; fv<cell.faces[f].vertex_ids.size(); fv++)
        {
          int fvi = cell.faces[f].vertex_ids[fv];
          nodal_boundary_numbers[fvi] =
            cell.faces[f].neighbor;
        }//for fv
      }//if bndry
    }//for face v's
  }//for local cell
  MPI_Barrier(MPI_COMM_WORLD);



  chi_log.Log(LOG_0) << "Communicating border cell information.";
  chi_log.Log(LOG_0) << "Serializing border cell information.";
  //================================================== Serialize local cells
  // The vectorized values will be as follows
  // - cell_glob_index
  // - cell_dof_start
  // - cell_type
  // - cell_mat_id
  // - cell_dof_count
  // - cell_face_count
  //
  // - dof 0 glob_index
  //     to
  // - dof N glob_index
  //
  // - face_0 dof_count
  // - face_0 dof 0 glob_index
  //     to
  // - face_0 dof fN glob_index
  //
  // - repeat all face info
  std::vector<int> border_cell_info;

  //============================================= Loop over set
  std::set<int>::iterator local_cell;
  for (local_cell  = local_border_cells.begin();
       local_cell != local_border_cells.end();
       local_cell++)
  {
    int local_cell_index = *local_cell;
    int cell_glob_index = grid->local_cell_glob_indices[local_cell_index];

    auto cell = &grid->local_cells[local_cell_index];
    DiffusionIPCellView* ip_view = ip_cell_views[local_cell_index];

    border_cell_info.push_back(cell_glob_index);         //cell_glob_index
    border_cell_info.push_back(ip_view->cell_dof_start); //cell_dof_start

    if (cell->Type() == chi_mesh::CellType::SLAB)
      border_cell_info.push_back(3);                         //cell_type
    if (cell->Type() == chi_mesh::CellType::POLYGON)
      border_cell_info.push_back(4);                         //cell_type
    if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
      border_cell_info.push_back(5);                         //cell_type

    border_cell_info.push_back(cell->material_id);      //cell_mat_id
    border_cell_info.push_back(cell->vertex_ids.size());//cell_dof_count
    border_cell_info.push_back(cell->faces.size());     //cell_face_count

    for (int v=0; v<cell->vertex_ids.size(); v++)
      border_cell_info.push_back(cell->vertex_ids[v]);//dof 0 to N

    for (int f=0; f<cell->faces.size(); f++)
    {
      int face_dof_count = cell->faces[f].vertex_ids.size();
      border_cell_info.push_back(face_dof_count);         //face dof_count
      for (int fv=0; fv<face_dof_count; fv++)
        border_cell_info.push_back(cell->faces[f].vertex_ids[fv]);
      //face dof 0 to fN
    }
  }//for local cell

  chi_log.Log(LOG_0) << "Broadcasting border cell information.";
  MPI_Barrier(MPI_COMM_WORLD);

  //================================================== Distribute border info
  std::vector<int> locI_info_size;
  std::vector<std::vector<int>> locI_border_cell_info;

  locI_info_size.resize(chi_mpi.process_count);
  locI_border_cell_info.resize(chi_mpi.process_count);

  //======================================== Collect sizes
  for (int locI=0; locI<chi_mpi.process_count; locI++)
  {
    if (locI == chi_mpi.location_id)
    {
      locI_info_size[locI] = border_cell_info.size();
    }
    MPI_Bcast(&locI_info_size[locI],1,MPI_INT,locI,MPI_COMM_WORLD);
  }

  //======================================== Collect info
  for (int locI=0; locI<chi_mpi.process_count; locI++)
  {
    if (locI == chi_mpi.location_id)
    {
      std::copy(border_cell_info.begin(),
                border_cell_info.end(),
                std::back_inserter(locI_border_cell_info[locI]));
    }
    else
      locI_border_cell_info[locI].resize(locI_info_size[locI]);

    MPI_Bcast(locI_border_cell_info[locI].data(),
              locI_info_size[locI],MPI_INT,locI,MPI_COMM_WORLD);
  }

  if (true)
    chi_log.Log(LOG_0) << "Deserializing border cell information.";
  MPI_Barrier(MPI_COMM_WORLD);

  //================================================== Deserialize border info
  // The vectorized values will be as follows
  // - cell_glob_index
  // - cell_dof_start
  // - cell_type
  // - cell_mat_id
  // - cell_dof_count
  // - cell_face_count
  //
  // - dof 0 glob_index
  //     to
  // - dof N glob_index
  //
  // - face_0 dof_count
  // - face_0 dof 0 glob_index
  //     to
  // - face_0 dof fN glob_index
  //
  // - repeat all face info
  ip_locI_bordercell_info.resize(chi_mpi.process_count);
  ip_locI_bordercells.resize(chi_mpi.process_count);
  ip_locI_borderfeviews.resize(chi_mpi.process_count);
  ip_locI_borderipviews.resize(chi_mpi.process_count);
  for (int locI=0; locI<chi_mpi.process_count; locI++)
  {
    int k=0;
    while (k<locI_info_size[locI])
    {
      DiffusionIPBorderCell* border_cell = new DiffusionIPBorderCell;
      border_cell->cell_glob_index = locI_border_cell_info[locI][k]; k++;
      border_cell->cell_dof_start  = locI_border_cell_info[locI][k]; k++;
      border_cell->cell_type       = locI_border_cell_info[locI][k]; k++;
      border_cell->cell_mat_id     = locI_border_cell_info[locI][k]; k++;
      border_cell->cell_dof_count  = locI_border_cell_info[locI][k]; k++;
      border_cell->cell_face_count = locI_border_cell_info[locI][k]; k++;

      int dof_count = border_cell->cell_dof_count;

      for (int v=0; v<dof_count; v++)
      {
        border_cell->v_indices.push_back(locI_border_cell_info[locI][k]);
        k++;
      }

      int face_count = border_cell->cell_face_count;

      for (int f=0; f<face_count; f++)
      {
        int face_dof_count = locI_border_cell_info[locI][k]; k++;
        border_cell->face_v_indices.emplace_back();
        for (int fv=0; fv<face_dof_count; fv++)
        {
          int vgi = locI_border_cell_info[locI][k]; k++;
          border_cell->face_v_indices[f].push_back(vgi);
        }
      }

      ip_locI_bordercell_info[locI].push_back(border_cell);
    }//while less than buffersize

    int locI_num_bordercells  = ip_locI_bordercell_info[locI].size();
    ip_locI_bordercells[locI].resize(locI_num_bordercells,nullptr);
    ip_locI_borderfeviews[locI].resize(locI_num_bordercells,nullptr);
    ip_locI_borderipviews[locI].resize(locI_num_bordercells,nullptr);
  }

  //================================================== Building off-diagonal
  //                                                   sparsity pattern
  dof_count = 0;
  for (auto& cell : grid->local_cells)
  {
    for (int v=0; v<cell.vertex_ids.size(); v++)
    {
      for (int f=0; f<cell.faces.size(); f++)
      {
        int neighbor = cell.faces[f].neighbor;
        bool is_bndry = grid->IsCellBndry(neighbor);
        bool is_local = cell.faces[f].IsNeighborLocal(grid);

        if ((not is_bndry) and (not is_local))
        {
          int adj_cell_partition_id = cell.faces[f].GetNeighborPartitionID(grid);
          auto adj_polyh_cell = (chi_mesh::Cell*)
            GetBorderCell(adj_cell_partition_id, neighbor);
          nodal_nnz_off_diag[dof_count] += adj_polyh_cell->vertex_ids.size();
        }
      }//for face
      dof_count++;
    }//for v
  }//for local cell
  MPI_Barrier(MPI_COMM_WORLD);


  chi_log.Log(LOG_0) << "Done creating DFEM sparsity pattern";

}