#include "chi_sweep.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>

#include "chi_SPDS.h"

#include <chi_mpi.h>
#include <chi_log.h>
#include <ChiTimer/chi_timer.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

extern double chi_global_timings[20];

//###################################################################
/**Populates the local sub-grid connection information for sweep orderings.*/
void chi_mesh::sweep_management::PopulateCellRelationships(
         chi_mesh::MeshContinuum *grid,
         chi_mesh::sweep_management::SPDS* sweep_order,
         std::vector<std::set<int>>& cell_dependencies,
         std::vector<std::set<int>>& cell_successors)
{
  size_t num_loc_cells = grid->local_cell_glob_indices.size();

  double tolerance = 1.0e-8;

  //============================================= Make directed connections
  for (size_t c=0; c<num_loc_cells; c++)
  {
    size_t cell_index = grid->local_cell_glob_indices[c];
    chi_mesh::Cell* cell = grid->cells[cell_index];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      auto slab_cell = (chi_mesh::CellSlab*)cell;

      int num_faces = 2;
      for (int f=0; f<num_faces; f++)
      {
        //======================================= Determine if the face
        //                                        is incident
        bool is_outgoing = false;
        double dot_normal = sweep_order->omega.Dot(slab_cell->face_normals[f]);
        if (dot_normal>(0.0+tolerance)) {is_outgoing = true;}

        //======================================= If outgoing determine if
        //                                        it is to a local cell
        if (is_outgoing)
        {
          int adj_cell_glob_index = slab_cell->edges[f];

          //================================if it is a cell and not bndry
          if (adj_cell_glob_index>=0)
          {
            auto adj_cell = grid->cells[adj_cell_glob_index];

            //========================= If it is not the current location
            if (adj_cell->partition_id == chi_mpi.location_id)
            {
              int adj_cell_local_index =
                grid->glob_cell_local_indices[adj_cell_glob_index];
//              boost::add_edge(c,adj_cell_local_index,G);

              cell_successors[c].insert(adj_cell->cell_local_id);
            }
            else
              sweep_order->AddLocalSuccessor(adj_cell->partition_id);

          }
        }
          //======================================= If not outgoing determine
          //                                        what it is dependent on
        else
        {
          int adj_cell_glob_index = slab_cell->edges[f];

          //================================if it is a cell and not bndry
          if (adj_cell_glob_index>=0)
          {
            auto adj_cell = grid->cells[adj_cell_glob_index];

            if (adj_cell->partition_id == chi_mpi.location_id)
              cell_dependencies[c].insert(adj_cell->cell_local_id);
            else
              sweep_order->AddLocalDependecy(adj_cell->partition_id);

          }
        }

      }//for face
    }
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    else if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)cell;

      for (int e=0; e<poly_cell->edges.size(); e++)
      {
        //======================================= Determine if the face
        //                                        is incident
        bool is_outgoing = false;
        double dot_normal = sweep_order->omega.Dot(poly_cell->edgenormals[e]);
        if (dot_normal>(0.0+tolerance)) {is_outgoing = true;}

        //======================================= If outgoing determine if
        //                                        it is to a local cell
        if (is_outgoing)
        {
          int adj_cell_glob_index = poly_cell->edges[e][2];

          //================================if it is a cell and not bndry
          if (adj_cell_glob_index>=0)
          {
            auto adj_cell = grid->cells[adj_cell_glob_index];

            //========================= If it is not the current location
            if (adj_cell->partition_id == chi_mpi.location_id)
            {
              int adj_cell_local_index =
                grid->glob_cell_local_indices[adj_cell_glob_index];
//              boost::add_edge(c,adj_cell_local_index,G);

              cell_successors[c].insert(adj_cell->cell_local_id);
            }
            else
              sweep_order->AddLocalSuccessor(adj_cell->partition_id);

          }
        }
          //======================================= If not outgoing determine
          //                                        what it is dependent on
        else
        {
          int adj_cell_glob_index = poly_cell->edges[e][2];

          //================================if it is a cell and not bndry
          if (adj_cell_glob_index>=0)
          {
            auto adj_cell = grid->cells[adj_cell_glob_index];

            if (adj_cell->partition_id == chi_mpi.location_id)
              cell_dependencies[c].insert(adj_cell->cell_local_id);
            else
              sweep_order->AddLocalDependecy(adj_cell->partition_id);

          }
        }

      }//for edge
    } //If polygon
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;

      for (int f=0; f<polyh_cell->faces.size(); f++)
      {
        //======================================= Determine if the face
        //                                        is incident
        bool is_outgoing = false;
        double dot_normal = sweep_order->omega.Dot(polyh_cell->faces[f]->geometric_normal);
        if (dot_normal>(0.0+tolerance)) {is_outgoing = true;}

        //======================================= If outgoing determine if
        //                                        it is to a local cell
        if (is_outgoing)
        {
          int adj_cell_glob_index = polyh_cell->faces[f]->face_indices[0];

          //================================if it is a cell and not bndry
          if (adj_cell_glob_index>=0)
          {
            auto adj_cell = grid->cells[adj_cell_glob_index];

            //========================= If it is in the current location
            if (adj_cell->partition_id == chi_mpi.location_id)
            {
              int adj_cell_local_index =
                grid->glob_cell_local_indices[adj_cell_glob_index];
//              boost::add_edge(c,adj_cell_local_index,G);

              cell_successors[c].insert(adj_cell->cell_local_id);
            }
            else
              sweep_order->AddLocalSuccessor(adj_cell->partition_id);
          }

        }
          //======================================= If not outgoing determine
          //                                        what it is dependent on
        else
        {
          int adj_cell_glob_index = polyh_cell->faces[f]->face_indices[0];

          //================================if it is a cell and not bndry
          if (adj_cell_glob_index>=0)
          {
            auto adj_cell = grid->cells[adj_cell_glob_index];

            if (adj_cell->partition_id == chi_mpi.location_id)
              cell_dependencies[c].insert(adj_cell->cell_local_id);
            else
              sweep_order->AddLocalDependecy(adj_cell->partition_id);

          }
        }

      }//for edge
    }

  }//for cell
}



