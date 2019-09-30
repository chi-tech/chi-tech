#include "../chi_mesh.h"
#include "chi_sweep.h"

#include "../MeshHandler/chi_meshhandler.h"
#include "../Cell/cell_slab.h"
#include "../Cell/cell_polyhedron.h"
#include "../Cell/cell_polygon.h"

#include "../MeshContinuum/chi_meshcontinuum.h"
#include "../VolumeMesher/chi_volumemesher.h"

#include "../../ChiGraph/chi_graph.h"

#include "chi_SPDS.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

#include "../../ChiConsole/chi_console.h"

extern ChiConsole chi_console;

#include <ChiTimer/chi_timer.h>

extern ChiTimer   chi_program_timer;

#include <ChiGraph/chi_directed_graph.h>

//###################################################################
/**Develops a sweep ordering for a given angle for locally owned
 * cells.*/
chi_mesh::SweepManagement::SPDS* chi_mesh::SweepManagement::
CreateSweepOrder(double polar, double azimuthal,
                 chi_mesh::MeshContinuum *vol_continuum,int number_of_groups,
                 bool allow_cycles)
{
  chi_mesh::SweepManagement::SPDS* sweep_order =
    new chi_mesh::SweepManagement::SPDS;
  sweep_order->polar     = polar;
  sweep_order->azimuthal = azimuthal;

  sweep_order->omega.x = sin(polar)*cos(azimuthal);
  sweep_order->omega.y = sin(polar)*sin(azimuthal);
  sweep_order->omega.z = cos(polar);

  sweep_order->grid    = vol_continuum;

  chi_mesh::MeshHandler* cur_handler = GetCurrentHandler();
  chi_mesh::VolumeMesher* mesher     = cur_handler->volume_mesher;

  double tolerance = 1.0e-8;

  //============================================= Compute direction vector
  chi_mesh::Vector omega;
  omega.x = sin(polar)*cos(azimuthal);
  omega.y = sin(polar)*sin(azimuthal);
  omega.z = cos(polar);
  if (chi_mpi.location_id == 0)
  {
    char buff[100];
    snprintf(buff, sizeof(buff), "Omega = %f,%f,%f\n",omega.x,omega.y,omega.z);

    chi_log.Log(LOG_0VERBOSE_1) << buff;
  }

  //============================================= Add all local item_id to graph
  // The cell added to the graph vertex will have
  // the same index as the cell's local index
  CHI_D_GRAPH G;
  size_t num_loc_cells = vol_continuum->local_cell_glob_indices.size();
  for (int c=0; c<num_loc_cells; c++)
  {
    boost::add_vertex(G);
  }



  //============================================= Make directed connections
  for (int c=0; c<num_loc_cells; c++)
  {
    size_t cell_index = vol_continuum->local_cell_glob_indices[c];
    chi_mesh::Cell* cell = vol_continuum->cells[cell_index];

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
        double dot_normal = omega.Dot(slab_cell->face_normals[f]);
        if (dot_normal>(0.0+tolerance)) {is_outgoing = true;}

        //======================================= If outgoing determine if
        //                                        it is to a local cell
        if (is_outgoing)
        {
          int adj_cell_glob_index = slab_cell->edges[f];

          //================================if it is a cell and not bndry
          if (adj_cell_glob_index>=0)
          {
            auto adj_cell = vol_continuum->cells[adj_cell_glob_index];

            //========================= If it is not the current location
            if (adj_cell->partition_id == chi_mpi.location_id)
            {
              int adj_cell_local_index =
                vol_continuum->glob_cell_local_indices[adj_cell_glob_index];
              boost::add_edge(c,adj_cell_local_index,G);
              //dependencies_per_cell[adj_cell_glob_index]->push_back(cell_index);
            }else
            {
              sweep_order->AddLocalSuccessor(adj_cell->partition_id);
            }

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
            auto adj_cell = vol_continuum->cells[adj_cell_glob_index];

            //========================= If it is not the current location
            if (adj_cell->partition_id != chi_mpi.location_id)
            {
              sweep_order->AddLocalDependecy(adj_cell->partition_id);
            }
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
        double dot_normal = omega.Dot(poly_cell->edgenormals[e]);
        if (dot_normal>(0.0+tolerance)) {is_outgoing = true;}

        //======================================= If outgoing determine if
        //                                        it is to a local cell
        if (is_outgoing)
        {
          int adj_cell_glob_index = poly_cell->edges[e][2];

          //================================if it is a cell and not bndry
          if (adj_cell_glob_index>=0)
          {
            auto adj_cell = vol_continuum->cells[adj_cell_glob_index];

            //========================= If it is not the current location
            if (adj_cell->partition_id == chi_mpi.location_id)
            {
              int adj_cell_local_index =
                vol_continuum->glob_cell_local_indices[adj_cell_glob_index];
              boost::add_edge(c,adj_cell_local_index,G);
              //dependencies_per_cell[adj_cell_glob_index]->push_back(cell_index);
            }else
            {
              sweep_order->AddLocalSuccessor(adj_cell->partition_id);
            }

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
            auto adj_cell = vol_continuum->cells[adj_cell_glob_index];

            //========================= If it is not the current location
            if (adj_cell->partition_id != chi_mpi.location_id)
            {
              sweep_order->AddLocalDependecy(adj_cell->partition_id);
            }
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
        double dot_normal = omega.Dot(polyh_cell->faces[f]->geometric_normal);
        if (dot_normal>(0.0+tolerance)) {is_outgoing = true;}

        //======================================= If outgoing determine if
        //                                        it is to a local cell
        if (is_outgoing)
        {
          int adj_cell_glob_index = polyh_cell->faces[f]->face_indices[0];

          //================================if it is a cell and not bndry
          if (adj_cell_glob_index>=0)
          {
            auto adj_cell = vol_continuum->cells[adj_cell_glob_index];

            //========================= If it is in the current location
            if (adj_cell->partition_id == chi_mpi.location_id)
            {
              int adj_cell_local_index =
                vol_continuum->glob_cell_local_indices[adj_cell_glob_index];
              boost::add_edge(c,adj_cell_local_index,G);
              //dependencies_per_cell[adj_cell_glob_index]->push_back(cell_index);
            } else
            {
              sweep_order->AddLocalSuccessor(adj_cell->partition_id);
            }
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
            auto adj_cell = vol_continuum->cells[adj_cell_glob_index];

            //========================= If it is not the current location
            if (adj_cell->partition_id != chi_mpi.location_id)
            {
              sweep_order->AddLocalDependecy(adj_cell->partition_id);
            }
          }
        }

      }//for edge
    }

  }//for cell



  //================================================== Generic topological
  //                                                   sorting
  typedef boost::graph_traits<CHI_D_GRAPH>::vertex_descriptor gVertex;

  boost::property_map<CHI_D_GRAPH, boost::vertex_index_t>::type
    index_map = get(boost::vertex_index, G);

  std::vector<gVertex> sorted_list;
  try{
    boost::topological_sort(G,std::back_inserter(sorted_list));
  }
  catch (const boost::bad_graph& exc)
  {
    chi_log.Log(LOG_ALLERROR)
    << "Cyclic local sweep ordering detected.";
    exit(EXIT_FAILURE);
  }


  //================================================== Generating sweep planes
  //The functionality of a sweep order allows for creating
  //multiple sweep planes but since this is a local sweep
  //we do not require to sort the sweep planes by the
  //their degrees and we can just use a single sweep plane.
  //Alternatively this code can be modified to allows this
  //but I can see no reason for this other than for
  //visualization.
  int i=0;
  chi_mesh::SweepManagement::SPLS* new_swp;
  for (std::vector<gVertex>::reverse_iterator ii=sorted_list.rbegin();
       ii!=sorted_list.rend();
       ii++)
  {
    if (i==0)
    {
      new_swp = new chi_mesh::SweepManagement::SPLS;
      sweep_order->spls=new_swp;
    }

    int cell_local_index = index_map[*ii];
    int cell_global_index =
      vol_continuum->local_cell_glob_indices[cell_local_index];
    new_swp->item_id.push_back(cell_global_index);
    i++;

    //if (i==1000) {i=0;}

  }



  G.clearing_graph();

  chi_log.Log(LOG_0VERBOSE_1)
  << chi_program_timer.GetTimeString()
  << " Communicating sweep dependencies.";

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create Task Dependency Graphs
  //All locations will send their dependencies
  //to the other locations
  int P = chi_mpi.process_count;
  std::vector<int> dependency_count_per_location(P,0);
  std::vector<std::vector<int>> global_dependencies(P, std::vector<int>());

  dependency_count_per_location[chi_mpi.location_id] =
    sweep_order->location_dependencies.size();

  //============================================= Broadcast location dep counts
  for (int locI=0; locI<P; locI++)
  {
    MPI_Bcast(&dependency_count_per_location[locI], //Buffer
              1, MPI_INT,                           //Count and type
              locI,                                 //Sending location
              MPI_COMM_WORLD);                      //Communicator
  }

  //============================================= Broadcast dependencies
  for (int locI=0; locI<P; locI++)
  {
    std::vector<int> locI_dependencies;
    if (locI == chi_mpi.location_id)
    {
      std::copy(sweep_order->location_dependencies.begin(),
                sweep_order->location_dependencies.end(),
                std::back_inserter(global_dependencies[locI]));
    } else
    {
      global_dependencies[locI].
        resize(dependency_count_per_location[locI],-1);
    }

    MPI_Bcast(global_dependencies[locI].data(), //Buffer
              dependency_count_per_location[locI],     //Count
              MPI_INT,                                 //Type
              locI,                                    //Sending location
              MPI_COMM_WORLD);                         //Communicator
  }

  //====================================== Filter dependencies for cycles
  if (allow_cycles)
  {
    for (int locI=0; locI<P; locI++)
    {
      for (int c=0; c<global_dependencies[locI].size(); c++)
      {
        int rlocI = global_dependencies[locI][c];
        if (rlocI<0)
          continue;
        for (int d=0; d<global_dependencies[rlocI].size(); d++)
        {
          if (global_dependencies[rlocI][d] == locI)
          {
            global_dependencies[locI][c]  = -1;
//            global_dependencies[rlocI][d] = -1;

            if (locI == chi_mpi.location_id)
            {
              std::vector<int>::iterator dependent_location =
                  std::find(sweep_order->location_dependencies.begin(),
                            sweep_order->location_dependencies.end(),
                            rlocI);
              sweep_order->location_dependencies.erase(dependent_location);
              sweep_order->delayed_location_dependencies.push_back(rlocI);
            }

            if (rlocI == chi_mpi.location_id)
            {
              sweep_order->delayed_location_successors.push_back(locI);
            }
          }//if cyclic dependency
        }//for rlocI dependency d
      }//for locI dependency c
    }//for locI
  }//if cycles allowed
  MPI_Barrier(MPI_COMM_WORLD);

  //====================================== Build task dependency graph
  CHI_D_GRAPH TDG;
  std::vector<int> glob_linear_sweep_order;
  std::vector<int> glob_sweep_order_rank;
  std::vector<int> glob_order_mapping;

  //================================= Add vertices to the graph
  for (int loc=0; loc<chi_mpi.process_count; loc++)
  {
    boost::add_vertex(TDG);
    glob_sweep_order_rank.push_back(-1);
    glob_order_mapping.push_back(loc);
  }

  //================================= Add dependencies
  for (int loc=0; loc<chi_mpi.process_count; loc++)
  {
    for (int dep=0; dep<global_dependencies[loc].size(); dep++)
    {
      if (global_dependencies[loc][dep]>=0)
        boost::add_edge(global_dependencies[loc][dep], loc, TDG);
    }
  }

  //================================= Generate topological sort
  boost::property_map<CHI_D_GRAPH, boost::vertex_index_t>::type
    glob_index_map = get(boost::vertex_index, TDG);

  std::vector<gVertex> glob_sorted_list;
  try {
    boost::topological_sort(TDG, std::back_inserter(glob_sorted_list));
  }
  catch (const boost::bad_graph& exc)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Cyclic global sweep ordering detected.";
    exit(EXIT_FAILURE);
  }

  //================================= Generate linear ordering
  // This step just puts the topological
  // sorting into a std::vector.
  for (std::vector<gVertex>::reverse_iterator ii=glob_sorted_list.rbegin();
       ii!=glob_sorted_list.rend();
       ii++)
  {
    int location_index = glob_index_map[*ii];
    glob_linear_sweep_order.push_back(location_index);
  }
  if (glob_linear_sweep_order.size() == 0)
  {
    chi_log.Log(LOG_0ERROR) << "Empty linear sweep ordering.";
    exit(EXIT_FAILURE);
  }


  //================================= Compute reorder mapping
  // This mapping allows us to punch in
  // the location id and find what its
  // id is in the TDG
  size_t num_ord = glob_linear_sweep_order.size();
  for (int k=0; k<num_ord; k++)
  {
    int loc = glob_linear_sweep_order[k];
    glob_order_mapping[loc] = k;
  }

  //================================= Determine sweep order ranks
  int abs_max_rank = 0;
  for (int k=0; k<num_ord; k++)
  {
    int loc = glob_linear_sweep_order[k];
    if (global_dependencies[loc].size() == 0)
    {
      glob_sweep_order_rank[k] = 0;
    }
    else
    {
      int max_rank = -1;
      for (int dep=0; dep<global_dependencies[loc].size(); dep++)
      {
        int dep_loc = global_dependencies[loc][dep];
        int dep_mapped_index = glob_order_mapping[dep_loc];

        if (glob_sweep_order_rank[dep_mapped_index] > max_rank)
        {
          max_rank = glob_sweep_order_rank[dep_mapped_index];
        }
      }
      glob_sweep_order_rank[k] = max_rank + 1;
      if ((max_rank + 1) > abs_max_rank)
      {
        abs_max_rank = max_rank + 1;
      }
    }
  }

  //================================= Generate TDG structure
  for (int r=0; r<=abs_max_rank; r++)
  {
    chi_mesh::SweepManagement::STDG* new_stdg =
      new chi_mesh::SweepManagement::STDG;
    sweep_order->global_sweep_planes.push_back(new_stdg);

    for (int k=0; k<num_ord; k++)
    {
      if (glob_sweep_order_rank[k] == r)
      {
        new_stdg->item_id.push_back(glob_linear_sweep_order[k]);
      }
    }
  }


  TDG.clear();





  MPI_Barrier(MPI_COMM_WORLD);






  return sweep_order;
}
