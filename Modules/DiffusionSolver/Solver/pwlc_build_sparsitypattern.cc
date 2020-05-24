#include "diffusion_solver.h"

#include <chi_log.h>
extern ChiLog& chi_log;

#include <chi_mpi.h>
extern ChiMPI& chi_mpi;

//###################################################################
/** Builds the sparsity pattern for PWLC spatial discretization.*/
//void chi_diffusion::Solver::PWLCBuildSparsityPattern()
//{
//  size_t num_local_cells = grid->local_cell_glob_indices.size();
//  for (int lc=0; lc < num_local_cells; lc++)
//  {
//    int glob_index = grid->local_cell_glob_indices[lc];
//
//    auto cell = grid->cells[glob_index];
//
//    for (int i=0; i < cell->vertex_ids.size(); i++)
//    {
//      int ir =  mesher->MapNode(cell->vertex_ids[i]);
//
//      if (ir<0)
//      {
//        chi_log.Log(LOG_ALLERROR)
//          << "ir Mapping error node " << cell->vertex_ids[i];
//        exit(EXIT_FAILURE);
//      }
//
//      //================================== Check if i is on boundary
//      for (int f=0; f < cell->faces.size(); f++)
//      {
//        if (cell->faces[f].neighbor < 0)
//        {
//          chi_mesh::CellFace& face = cell->faces[f];
//          size_t num_face_verts = face.vertex_ids.size();
//          for (int fv=0; fv<num_face_verts; fv++)
//          {
//            int v0_index =
//              mesher->MapNode(face.vertex_ids[fv]);
//
//            if (v0_index<0)
//            {
//              chi_log.Log(LOG_ALLERROR)
//                << "v0 Mapping error node " << face.vertex_ids[fv];
//              exit(EXIT_FAILURE);
//            }
//
//            if (ir == v0_index)
//            {
//              //================= Processing boundary
//              int boundary_type =
//                boundaries[abs(face.neighbor)-1]->type;
//              if (boundary_type == DIFFUSION_DIRICHLET)
//              {
//                nodal_boundary_numbers[ir]= face.neighbor;
//              }
//              break;
//            } //if ir part of face
//          }
//        }
//      }//for f
//
//      //======================================= Set nodal connections
//      std::vector<int>& node_links = nodal_connections[ir];
//      for (int j=0; j < cell->vertex_ids.size(); j++)
//      {
//        int jr = mesher->MapNode(cell->vertex_ids[j]);
//
//        //====================== Check for duplicates
//        bool already_there = false;
//        for (int k=0; k<node_links.size(); k++)
//        {
//          if (node_links[k] == jr)
//          {already_there = true; break;}
//        }
//        if (!already_there)
//        {
//          node_links.push_back(jr);
//          if ((jr>=local_rows_from) && (jr<=local_rows_to))
//          {
//            nodal_nnz_in_diag[ir]+=1;
//          } else
//          {
//            nodal_nnz_off_diag[ir]+=1;
//          }
//        }
//      }//for j
//    }//for i
//  }
//
//  MPI_Barrier(MPI_COMM_WORLD);
//}