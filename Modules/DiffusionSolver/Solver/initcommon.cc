#include "diffusion_solver.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"
#include "ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h"
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include "../Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "../Boundaries/chi_diffusion_bndry_reflecting.h"

#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>
#include <ChiMesh/Cell/cell_newbase.h>

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;




//###################################################################
/**Initialization of common to all solver types.*/
void chi_diffusion::Solver::InitializeCommonItems()
{
  chi_mesh::Region*  region = this->regions.back();

  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherExtruder))
  {
    for (int b=0; b<(region->boundaries.size()-2); b++)
    {
      chi_diffusion::Boundary* new_bndry =
        new chi_diffusion::BoundaryDirichlet;
      this->boundaries.push_back(new_bndry);
      chi_log.Log(LOG_0VERBOSE_1)
      << "Dirichlet boundary added (index " << boundaries.size()-1 <<  ").";
    }
    chi_diffusion::Boundary* new_bndry;

    new_bndry = new chi_diffusion::BoundaryReflecting;
    this->boundaries.push_back(new_bndry);
    chi_log.Log(LOG_0VERBOSE_1)
      << "Reflecting boundary added (index " << boundaries.size()-1 <<  ").";

    new_bndry = new chi_diffusion::BoundaryReflecting;
    this->boundaries.push_back(new_bndry);
    chi_log.Log(LOG_0VERBOSE_1)
      << "Reflecting boundary added (index " << boundaries.size()-1 <<  ").";
  }
  else if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
  {
    chi_diffusion::Boundary* new_bndry =
      new chi_diffusion::BoundaryDirichlet;
    this->boundaries.push_back(new_bndry);
    chi_log.Log(LOG_0VERBOSE_1)
      << "Dirichlet boundary added (index " << boundaries.size()-1 <<  ").";

    new_bndry =
      new chi_diffusion::BoundaryDirichlet;
    this->boundaries.push_back(new_bndry);
    chi_log.Log(LOG_0VERBOSE_1)
      << "Dirichlet boundary added (index " << boundaries.size()-1 <<  ").";
  }
  else
  {
    for (int b=0; b<region->boundaries.size(); b++)
    {
      chi_diffusion::Boundary* new_bndry =
        new chi_diffusion::BoundaryDirichlet;
      this->boundaries.push_back(new_bndry);
      chi_log.Log(LOG_0VERBOSE_1)
        << "Dirichlet boundary added (index " << boundaries.size()-1 <<  ").";
    }
  }

  common_items_initialized = true;
}


//###################################################################
/**Distributes border cell information.*/
void chi_diffusion::Solver::PWLDBuildSparsityPattern()
{
  int num_loc_cells = grid->local_cell_glob_indices.size();
  int dof_count = 0;
  std::set<int> local_border_cells;
  for (int lc=0; lc<num_loc_cells; lc++)
  {
    int cell_glob_index = grid->local_cell_glob_indices[lc];
    auto cell = grid->cells[cell_glob_index];

//    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
//    if (cell->Type() == chi_mesh::CellType::SLAB)
//    {
//      auto slab_cell = (chi_mesh::CellSlab*)cell;
//
//      DiffusionIPCellView* ip_view = new DiffusionIPCellView;
//      ip_view->cell_dof_start = dof_count + pwld_local_dof_start;
//      pwld_cell_dof_array_address.push_back(dof_count);
//      ip_cell_views.push_back(ip_view);
//
//      nodal_nnz_in_diag[dof_count]   = 4;
//      nodal_nnz_in_diag[dof_count+1] = 4;
//
//      if (slab_cell->edges[0]>=0)
//      {
//        if (grid->IsCellLocal(slab_cell->edges[0]))
//          nodal_nnz_in_diag[dof_count]  += 2;
//        else
//        {
////          nodal_nnz_off_diag[dof_count] += 2;
//          local_border_cells.insert(lc);
//        }
//      }
//      else
//        nodal_boundary_numbers[slab_cell->v_indices[0]] = slab_cell->edges[0];
//
//      if (slab_cell->edges[1]>=0)
//      {
//        if (grid->IsCellLocal(slab_cell->edges[1]))
//          nodal_nnz_in_diag[dof_count+1]  += 2;
//        else
//        {
////          nodal_nnz_off_diag[dof_count+1] += 2;
//          local_border_cells.insert(lc);
//        }
//      }
//      else
//        nodal_boundary_numbers[slab_cell->v_indices[1]] = slab_cell->edges[1];
//
//      dof_count += 2;
//    }
//
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)cell;

      DiffusionIPCellView* ip_view = new DiffusionIPCellView;
      ip_view->cell_dof_start = dof_count + pwld_local_dof_start;
      pwld_cell_dof_array_address.push_back(dof_count);
      ip_cell_views.push_back(ip_view);

      for (int v=0; v<poly_cell->v_indices.size(); v++)
      {
        nodal_nnz_in_diag[dof_count] = poly_cell->v_indices.size();

        for (int e=0; e<poly_cell->edges.size(); e++)
        {
          if (poly_cell->edges[e][2]>=0) //Not boundary
          {
            if (grid->IsCellLocal(poly_cell->edges[e][2]))
            {
              int adj_cell_glob_index = poly_cell->edges[e][2];
              auto adj_cell =
                (chi_mesh::CellPolygon*)grid->cells[adj_cell_glob_index];
              nodal_nnz_in_diag[dof_count] += adj_cell->v_indices.size();
            }
            else
            {
              //Since we have no information about the non-local cell,
              //we can make a good assumption that it has the same amount
              //of dofs than does the current cell.
//              nodal_nnz_off_diag[dof_count] += 2*poly_cell->v_indices.size();
              local_border_cells.insert(lc);
            }
          }//if not bndry
        }//for edge
        dof_count++;
      }//for vi

      //==================================== Boundary numbers
      for (int e=0; e<poly_cell->edges.size(); e++)
      {
        if (poly_cell->edges[e][2]<0)
        {
          nodal_boundary_numbers[poly_cell->edges[e][0]]=
            poly_cell->edges[e][2];
          nodal_boundary_numbers[poly_cell->edges[e][1]]=
            poly_cell->edges[e][2];
        }
      }

    }//polygon

//    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
//    if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
//    {
//      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
//
//      DiffusionIPCellView* ip_view = new DiffusionIPCellView;
//      ip_view->cell_dof_start = dof_count + pwld_local_dof_start;
//      pwld_cell_dof_array_address.push_back(dof_count);
//      ip_cell_views.push_back(ip_view);
//
//      for (int v=0; v<polyh_cell->v_indices.size(); v++)
//      {
//        nodal_nnz_in_diag[dof_count] = polyh_cell->v_indices.size();
//
//        for (int f=0; f<polyh_cell->faces.size(); f++)
//        {
//          if (polyh_cell->faces[f]->face_indices[NEIGHBOR]>=0) //Not bndry
//          {
//            bool is_local =
//              grid->IsCellLocal(polyh_cell->faces[f]->face_indices[NEIGHBOR]);
//
//            if (is_local)
//            {
//              int adj_cell_glob_index =
//                polyh_cell->faces[f]->face_indices[NEIGHBOR];
//              auto adj_cell =
//                (chi_mesh::CellPolyhedron*)grid->cells[adj_cell_glob_index];
//              nodal_nnz_in_diag[dof_count] += adj_cell->v_indices.size();
//            }
//            else
//            {
//              //Since we have no information about the non-local cell,
//              //we can make a good assumption that it has the same amount
//              //of dofs than does the current cell.
////              nodal_nnz_off_diag[dof_count] += polyh_cell->v_indices.size();
//              local_border_cells.insert(lc);
//            }
//          }
//        }
//        dof_count++;
//      }
//
//      //==================================== Boundary numbers
//      for (int f=0; f<polyh_cell->faces.size(); f++)
//      {
//        if (polyh_cell->faces[f]->face_indices[NEIGHBOR]<0)
//        {
//          for (int fv=0; fv<polyh_cell->faces[f]->v_indices.size(); fv++)
//          {
//            int fvi = polyh_cell->faces[f]->v_indices[fv];
//            nodal_boundary_numbers[fvi] =
//              polyh_cell->faces[f]->face_indices[NEIGHBOR];
//          }//for fv
//        }//if bndry
//      }//for face v's
//    }//if polyhedron

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CELL_NEWBASE
    if (cell->Type() == chi_mesh::CellType::CELL_NEWBASE)
    {
      auto cell_base = (chi_mesh::CellBase*)cell;

      DiffusionIPCellView* ip_view = new DiffusionIPCellView;
      ip_view->cell_dof_start = dof_count + pwld_local_dof_start;
      pwld_cell_dof_array_address.push_back(dof_count);
      ip_cell_views.push_back(ip_view);

      for (int v=0; v<cell_base->vertex_ids.size(); v++)
      {
        nodal_nnz_in_diag[dof_count] = cell_base->vertex_ids.size();

        for (int f=0; f<cell_base->faces.size(); f++)
        {
          if (cell_base->faces[f].neighbor >= 0) //Not bndry
          {
            bool is_local = grid->IsCellLocal(cell_base->faces[f].neighbor);

            if (is_local)
            {
              int adj_cell_glob_index = cell_base->faces[f].neighbor;
              auto adj_cell =
                (chi_mesh::CellBase*)grid->cells[adj_cell_glob_index];
              nodal_nnz_in_diag[dof_count] += adj_cell->vertex_ids.size();
            }
            else
            {
              local_border_cells.insert(lc);
            }
          }
        }
        dof_count++;
      }

      //==================================== Boundary numbers
      for (int f=0; f<cell_base->faces.size(); f++)
      {
        if (cell_base->faces[f].neighbor < 0)
        {
          for (int fv=0; fv<cell_base->faces[f].vertex_ids.size(); fv++)
          {
            int fvi = cell_base->faces[f].vertex_ids[fv];
            nodal_boundary_numbers[fvi] =
              cell_base->faces[f].neighbor;
          }//for fv
        }//if bndry
      }//for face v's
    }//if polyhedron

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

    auto cell = grid->cells[cell_glob_index];
    DiffusionIPCellView* ip_view = ip_cell_views[local_cell_index];

    border_cell_info.push_back(cell_glob_index);         //cell_glob_index
    border_cell_info.push_back(ip_view->cell_dof_start); //cell_dof_start

//    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
//    if (cell->Type() == chi_mesh::CellType::SLAB)
//    {
//      auto slab_cell = (chi_mesh::CellSlab*)cell;
//      border_cell_info.push_back(0);                     //cell_type
//      border_cell_info.push_back(slab_cell->material_id);//cell_mat_id
//      border_cell_info.push_back(2);                     //cell_dof_count
//      border_cell_info.push_back(2);                     //cell_face_count
//
//      for (int v=0; v<2; v++)
//        border_cell_info.push_back(slab_cell->v_indices[v]); //dof 0 to N
//
//      for (int f=0; f<2; f++)
//      {
//        border_cell_info.push_back(1);                   //face dof_count
//        border_cell_info.push_back(slab_cell->v_indices[f]); //face dof 0 to fN
//      }
//    }//slab
//
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)cell;
      border_cell_info.push_back(1);                     //cell_type
      border_cell_info.push_back(poly_cell->material_id);//cell_mat_id
      border_cell_info.push_back(poly_cell->v_indices.size());//cell_dof_count
      border_cell_info.push_back(poly_cell->edges.size());//cell_face_count

      for (int v=0; v<poly_cell->v_indices.size(); v++)
        border_cell_info.push_back(poly_cell->v_indices[v]);//dof 0 to N

      for (int f=0; f<poly_cell->edges.size(); f++)
      {
        border_cell_info.push_back(2);                       //face dof_count
        border_cell_info.push_back(poly_cell->edges[f][0]); //face dof 0 to fN
        border_cell_info.push_back(poly_cell->edges[f][1]); //face dof 0 to fN
      }
    }//polygon
//
//    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
//    if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
//    {
//      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
//      border_cell_info.push_back(2);                     //cell_type
//      border_cell_info.push_back(polyh_cell->material_id);//cell_mat_id
//      border_cell_info.push_back(polyh_cell->v_indices.size());//cell_dof_count
//      border_cell_info.push_back(polyh_cell->faces.size());//cell_face_count
//
//      for (int v=0; v<polyh_cell->v_indices.size(); v++)
//        border_cell_info.push_back(polyh_cell->v_indices[v]);//dof 0 to N
//
//      for (int f=0; f<polyh_cell->faces.size(); f++)
//      {
//        int face_dof_count = polyh_cell->faces[f]->v_indices.size();
//        border_cell_info.push_back(face_dof_count);         //face dof_count
//        for (int fv=0; fv<face_dof_count; fv++)
//          border_cell_info.push_back(polyh_cell->faces[f]->v_indices[fv]);
//        //face dof 0 to fN
//      }
//    }//polyhedron

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CELL_NEWBASE
    if (cell->Type() == chi_mesh::CellType::CELL_NEWBASE)
    {
      auto cell_base = (chi_mesh::CellBase*)cell;
      if (cell_base->Type2() == chi_mesh::CellType::POLYGONV2)
        border_cell_info.push_back(4);                         //cell_type
      if (cell_base->Type2() == chi_mesh::CellType::POLYHEDRONV2)
        border_cell_info.push_back(5);                         //cell_type

      border_cell_info.push_back(cell_base->material_id);      //cell_mat_id
      border_cell_info.push_back(cell_base->vertex_ids.size());//cell_dof_count
      border_cell_info.push_back(cell_base->faces.size());     //cell_face_count

      for (int v=0; v<cell_base->vertex_ids.size(); v++)
        border_cell_info.push_back(cell_base->vertex_ids[v]);//dof 0 to N

      for (int f=0; f<cell_base->faces.size(); f++)
      {
        int face_dof_count = cell_base->faces[f].vertex_ids.size();
        border_cell_info.push_back(face_dof_count);         //face dof_count
        for (int fv=0; fv<face_dof_count; fv++)
          border_cell_info.push_back(cell_base->faces[f].vertex_ids[fv]);
        //face dof 0 to fN
      }
    }//new cell base
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
  for (int lc=0; lc<num_loc_cells; lc++)
  {
    int cell_glob_index = grid->local_cell_glob_indices[lc];
    auto cell = grid->cells[cell_glob_index];

//    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
//    if (cell->Type() == chi_mesh::CellType::SLAB)
//    {
//      auto slab_cell = (chi_mesh::CellSlab*)cell;
//
//      nodal_nnz_in_diag[dof_count]   += 2;
//      nodal_nnz_in_diag[dof_count+1] += 2;
//      dof_count += 2;
//    }
//
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)cell;

      for (int v=0; v<poly_cell->v_indices.size(); v++)
      {
        for (int e=0; e<poly_cell->edges.size(); e++)
        {
          int neighbor = poly_cell->edges[e][EDGE_NEIGHBOR];
          bool is_bndry = grid->IsCellBndry(neighbor);
          bool is_local = grid->IsCellLocal(neighbor);

          if ((not is_bndry) and (not is_local))
          {
            auto adj_cell = grid->cells[neighbor];
            auto adj_poly_cell = (chi_mesh::CellPolygon*)
              GetBorderCell(adj_cell->partition_id,
                            neighbor);
            nodal_nnz_off_diag[dof_count] += adj_poly_cell->v_indices.size();
          }
        }//for edge
        dof_count++;
      }//for vi

    }//polygon
//
//    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
//    if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
//    {
//      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
//
//      for (int v=0; v<polyh_cell->v_indices.size(); v++)
//      {
//        for (int f=0; f<polyh_cell->faces.size(); f++)
//        {
//          int neighbor = polyh_cell->faces[f]->face_indices[NEIGHBOR];
//          bool is_bndry = grid->IsCellBndry(neighbor);
//          bool is_local = grid->IsCellLocal(neighbor);
//
//          if ((not is_bndry) and (not is_local))
//          {
//            auto adj_cell = grid->cells[neighbor];
//            auto adj_polyh_cell = (chi_mesh::CellPolyhedron*)
//                                  GetBorderCell(adj_cell->partition_id,
//                                                neighbor);
//            nodal_nnz_off_diag[dof_count] += adj_polyh_cell->v_indices.size();
//          }
//        }//for face
//        dof_count++;
//      }
//    }//if polyhedron

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CELL_NEWBASE
    if (cell->Type() == chi_mesh::CellType::CELL_NEWBASE)
    {
      auto cell_base = (chi_mesh::CellBase*)cell;

      for (int v=0; v<cell_base->vertex_ids.size(); v++)
      {
        for (int f=0; f<cell_base->faces.size(); f++)
        {
          int neighbor = cell_base->faces[f].neighbor;
          bool is_bndry = grid->IsCellBndry(neighbor);
          bool is_local = grid->IsCellLocal(neighbor);

          if ((not is_bndry) and (not is_local))
          {
            auto adj_cell = grid->cells[neighbor];
            auto adj_polyh_cell = (chi_mesh::CellBase*)
              GetBorderCell(adj_cell->partition_id,
                            neighbor);
            nodal_nnz_off_diag[dof_count] += adj_polyh_cell->vertex_ids.size();
          }
        }//for face
        dof_count++;
      }
    }//if polyhedron

  }//for local cell
  MPI_Barrier(MPI_COMM_WORLD);




}