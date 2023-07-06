#include "SPDS.h"

#include "graphs/chi_directed_graph.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"
#include "console/chi_console.h"
#include "utils/chi_timer.h"

#include <algorithm>

// ###################################################################
/** Given a location J index, maps to a predecessor location.*/
int chi_mesh::sweep_management::SPDS::MapLocJToPrelocI(int locJ) const
{
  for (int i = 0; i < location_dependencies_.size(); i++)
  {
    if (location_dependencies_[i] == locJ) { return i; }
  }

  for (int i = 0; i < delayed_location_dependencies_.size(); i++)
  {
    if (delayed_location_dependencies_[i] == locJ) { return -(i + 1); }
  }

  Chi::log.LogAllError()
    << "SPDS Invalid mapping encountered in MapLocJToPrelocI.";
  Chi::Exit(EXIT_FAILURE);
  return 0;
}

// ###################################################################
/** Given a location J index, maps to a dependent location.*/
int chi_mesh::sweep_management::SPDS::MapLocJToDeplocI(int locJ) const
{
  for (int i = 0; i < location_successors_.size(); i++)
  {
    if (location_successors_[i] == locJ) { return i; }
  }

  Chi::log.LogAllError()
    << "SPDS Invalid mapping encountered in MapLocJToDeplocI.";
  Chi::Exit(EXIT_FAILURE);
  return 0;
}

// ###################################################################
/**Populates cell relationships*/
void chi_mesh::sweep_management::SPDS::PopulateCellRelationships(
  const chi_mesh::Vector3& omega,
  std::set<int>& location_dependencies,
  std::set<int>& location_successors,
  std::vector<std::set<std::pair<int, double>>>& cell_successors)
{
  constexpr double tolerance = 1.0e-16;

  constexpr auto FOPARALLEL = FaceOrientation::PARALLEL;
  constexpr auto FOINCOMING = FaceOrientation::INCOMING;
  constexpr auto FOOUTGOING = FaceOrientation::OUTGOING;

  cell_face_orientations_.assign(grid_.local_cells.size(), {});
  for (auto& cell : grid_.local_cells)
    cell_face_orientations_[cell.local_id_].assign(cell.faces_.size(),
                                                   FOPARALLEL);

  for (auto& cell : grid_.local_cells)
  {
    size_t f = 0;
    for (auto& face : cell.faces_)
    {
      //======================================= Determine if the face
      //                                        is incident
      FaceOrientation orientation = FOPARALLEL;
      const double mu = omega.Dot(face.normal_);

      bool owns_face = true;
      if (face.has_neighbor_ and cell.global_id_ > face.neighbor_id_ and
          grid_.IsCellLocal(face.neighbor_id_))
        owns_face = false;

      if (owns_face)
      {
        // clang-format off
        if (mu > tolerance) orientation = FOOUTGOING;
        else if (mu < tolerance) orientation = FOINCOMING;

        cell_face_orientations_[cell.local_id_][f] = orientation;

        if (face.has_neighbor_ and grid_.IsCellLocal(face.neighbor_id_))
        {
          const auto& adj_cell = grid_.cells[face.neighbor_id_];
          const auto ass_face = face.GetNeighborAssociatedFace(grid_);
          auto& adj_face_ori =
            cell_face_orientations_[adj_cell.local_id_][ass_face];

          switch (orientation)
          {
            case FOPARALLEL: adj_face_ori = FOPARALLEL; break;
            case FOINCOMING: adj_face_ori = FOOUTGOING; break;
            case FOOUTGOING: adj_face_ori = FOINCOMING; break;
          }
        }
        // clang-format on
      } // if face owned
      else if (face.has_neighbor_ and not grid_.IsCellLocal(face.neighbor_id_))
      {
        const auto& adj_cell = grid_.cells[face.neighbor_id_];
        const auto ass_face = face.GetNeighborAssociatedFace(grid_);
        const auto& adj_face = adj_cell.faces_[ass_face];

        auto& cur_face_ori = cell_face_orientations_[cell.local_id_][f];

        const double adj_mu = omega.Dot(adj_face.normal_);
        if (adj_mu > tolerance) orientation = FOOUTGOING;
        else if (adj_mu < tolerance)
          orientation = FOINCOMING;

        switch (orientation)
        {
          case FOPARALLEL:
            cur_face_ori = FOPARALLEL;
            break;
          case FOINCOMING:
            cur_face_ori = FOOUTGOING;
            break;
          case FOOUTGOING:
            cur_face_ori = FOINCOMING;
            break;
        }
      } // if not face owned locally at all

      ++f;
    } // for face
  }

  //============================================= Make directed connections
  for (auto& cell : grid_.local_cells)
  {
    const uint64_t c = cell.local_id_;
    size_t f = 0;
    for (auto& face : cell.faces_)
    {
      const double mu = omega.Dot(face.normal_);
      //======================================= If outgoing determine if
      //                                        it is to a local cell
      if (cell_face_orientations_[cell.local_id_][f] == FOOUTGOING)
      {
        //================================ If it is a cell and not bndry
        if (face.has_neighbor_)
        {
          //========================= If it is in the current location
          if (face.IsNeighborLocal(grid_))
          {
            double weight = mu * face.ComputeFaceArea(grid_);
            cell_successors[c].insert(
              std::make_pair(face.GetNeighborLocalID(grid_), weight));
          }
          else
            location_successors.insert(face.GetNeighborPartitionID(grid_));
        }
      }
      //======================================= If not outgoing determine
      //                                        what it is dependent on
      else
      {
        //================================if it is a cell and not bndry
        if (face.has_neighbor_ and not face.IsNeighborLocal(grid_))
          location_dependencies.insert(face.GetNeighborPartitionID(grid_));
      }
      ++f;
    } // for face
  }   // for cell
}

// ###################################################################
/**Builds the task dependency graph.*/
void chi_mesh::sweep_management::SPDS::BuildTaskDependencyGraph(
  const std::vector<std::vector<int>>& global_dependencies,
  bool cycle_allowance_flag)
{

  std::vector<std::pair<int, int>> edges_to_remove;
  std::vector<int> raw_edges_to_remove;
  chi::DirectedGraph TDG;

  //============================================= Build graph on home location
  if (Chi::mpi.location_id == 0)
  {
    Chi::log.Log0Verbose1() << Chi::program_timer.GetTimeString()
                            << " Building Task Dependency Graphs.";

    //====================================== Add vertices to the graph
    for (int loc = 0; loc < Chi::mpi.process_count; loc++)
      TDG.AddVertex();

    //====================================== Add dependencies
    for (int loc = 0; loc < Chi::mpi.process_count; loc++)
      for (int dep = 0; dep < global_dependencies[loc].size(); dep++)
        TDG.AddEdge(global_dependencies[loc][dep], loc);

    //====================================== Remove cyclic dependencies
    if (cycle_allowance_flag)
    {
      Chi::log.Log0Verbose1() << Chi::program_timer.GetTimeString()
                              << " Removing intra-cellset cycles.";
      auto edges_to_remove_temp = TDG.RemoveCyclicDependencies();
      for (const auto& [v0, v1] : edges_to_remove_temp)
        edges_to_remove.emplace_back(v0, v1);
    }

    //====================================== Serialize edges to be removed
    raw_edges_to_remove.resize(edges_to_remove.size() * 2, 0);
    int i = 0;
    for (const auto& edge : edges_to_remove)
    {
      raw_edges_to_remove[i++] = edge.first;
      raw_edges_to_remove[i++] = edge.second;
    }
  } // if home

  //============================================= Broadcast edge buffer size
  int edge_buffer_size = 0;

  if (Chi::mpi.location_id == 0)
    edge_buffer_size = static_cast<int>(raw_edges_to_remove.size());

  MPI_Bcast(&edge_buffer_size, // Buffer
            1,
            MPI_INT,        // Count and datatype
            0,              // Root location
            Chi::mpi.comm); // Communicator

  //============================================= Broadcast edges
  if (Chi::mpi.location_id != 0)
    raw_edges_to_remove.resize(edge_buffer_size, -1);

  MPI_Bcast(raw_edges_to_remove.data(), // Buffer
            edge_buffer_size,
            MPI_INT,        // Count and datatype
            0,              // Root location
            Chi::mpi.comm); // Communicator

  //============================================= De-serialize edges
  if (Chi::mpi.location_id != 0)
  {
    edges_to_remove.resize(edge_buffer_size / 2, std::pair<int, int>(0, 0));
    int i = 0;
    for (auto& edge : edges_to_remove)
    {
      edge.first = raw_edges_to_remove[i++];
      edge.second = raw_edges_to_remove[i++];
    }
  }

  //============================================= Remove edges
  for (auto& edge_to_remove : edges_to_remove)
  {
    int rlocI = edge_to_remove.first;
    int locI = edge_to_remove.second;

    if (Chi::mpi.location_id == 0) TDG.RemoveEdge(rlocI, locI);

    if (locI == Chi::mpi.location_id)
    {
      auto dependent_location = std::find(
        location_dependencies_.begin(), location_dependencies_.end(), rlocI);
      location_dependencies_.erase(dependent_location);
      delayed_location_dependencies_.push_back(rlocI);
    }

    if (rlocI == Chi::mpi.location_id)
      delayed_location_successors_.push_back(locI);
  }

  //============================================= Generate topological sort
  std::vector<int> glob_linear_sweep_order;
  if (Chi::mpi.location_id == 0)
  {
    Chi::log.LogAllVerbose2() << Chi::program_timer.GetTimeString()
                              << "   - Generating topological sort.";
    auto so_temp = TDG.GenerateTopologicalSort();
    for (auto v : so_temp)
      glob_linear_sweep_order.emplace_back(v);

    if (glob_linear_sweep_order.empty())
    {
      Chi::log.LogAllError()
        << "Topological sorting for global sweep-ordering failed. "
        << "Cyclic dependencies detected. Cycles need to be allowed"
        << " by calling application.";
      Chi::Exit(EXIT_FAILURE);
    }
  }

  //============================================= Broadcasting topsort size
  int topsort_buffer_size = 0;

  if (Chi::mpi.location_id == 0)
    topsort_buffer_size = glob_linear_sweep_order.size();

  MPI_Bcast(&topsort_buffer_size, // Buffer
            1,
            MPI_INT,        // Count and datatype
            0,              // Root location
            Chi::mpi.comm); // Communicator

  //============================================= Broadcast topological sort
  if (Chi::mpi.location_id != 0)
    glob_linear_sweep_order.resize(topsort_buffer_size, -1);

  MPI_Bcast(glob_linear_sweep_order.data(), // Buffer
            topsort_buffer_size,
            MPI_INT,        // Count and datatype
            0,              // Root location
            Chi::mpi.comm); // Communicator

  //============================================= Compute reorder mapping
  // This mapping allows us to punch in
  // the location id and find what its
  // id is in the TDG
  std::vector<int> glob_order_mapping(Chi::mpi.process_count, -1);

  for (int k = 0; k < Chi::mpi.process_count; k++)
  {
    int loc = glob_linear_sweep_order[k];
    glob_order_mapping[loc] = k;
  }

  //============================================= Determine sweep order ranks
  Chi::log.Log0Verbose1() << Chi::program_timer.GetTimeString()
                          << " Determining sweep order ranks.";

  std::vector<int> glob_sweep_order_rank(Chi::mpi.process_count, -1);

  int abs_max_rank = 0;
  for (int k = 0; k < Chi::mpi.process_count; k++)
  {
    int loc = glob_linear_sweep_order[k];
    if (global_dependencies[loc].empty()) glob_sweep_order_rank[k] = 0;
    else
    {
      int max_rank = -1;
      for (auto dep_loc : global_dependencies[loc])
      {
        if (dep_loc < 0) continue;
        int dep_mapped_index = glob_order_mapping[dep_loc];

        if (glob_sweep_order_rank[dep_mapped_index] > max_rank)
          max_rank = glob_sweep_order_rank[dep_mapped_index];
      }
      glob_sweep_order_rank[k] = max_rank + 1;
      if ((max_rank + 1) > abs_max_rank) abs_max_rank = max_rank + 1;
    }
  }

  //============================================= Generate TDG structure
  Chi::log.Log0Verbose1() << Chi::program_timer.GetTimeString()
                          << " Generating TDG structure.";
  for (int r = 0; r <= abs_max_rank; r++)
  {
    chi_mesh::sweep_management::STDG new_stdg;

    for (int k = 0; k < Chi::mpi.process_count; k++)
    {
      if (glob_sweep_order_rank[k] == r)
        new_stdg.item_id.push_back(glob_linear_sweep_order[k]);
    }
    global_sweep_planes_.push_back(new_stdg);
  }
}

// ###################################################################
void chi_mesh::sweep_management::SPDS::PrintedGhostedGraph() const
{
  constexpr double tolerance = 1.0e-16;

  for (int p = 0; p < Chi::mpi.process_count; ++p)
  {
    if (p == Chi::mpi.location_id)
    {
      std::cout << "// Printing directed graph for location " << p << ":\n";
      std::cout << "digraph DG {\n";
      std::cout << "  splines=\"FALSE\";\n";
      std::cout << "  rankdir=\"LR\";\n\n";
      std::cout << "  /* Vertices */\n";

      for (const auto& cell : grid_.local_cells)
      {
        std::cout << "  " << cell.global_id_ << " [shape=\"circle\"]\n";

        for (const auto& face : cell.faces_)
        {
          if (face.has_neighbor_ and (not grid_.IsCellLocal(face.neighbor_id_)))
            std::cout << "  " << face.neighbor_id_
                      << " [shape=\"circle\", style=filled fillcolor=red "
                         "fontcolor=white] //proc="
                      << grid_.cells[face.neighbor_id_].partition_id_ << "\n";
        }
      }

      std::cout << "\n"
                << "  /* Edges */\n";
      for (const auto& cell : grid_.local_cells)
      {
        for (const auto& face : cell.faces_)
        {
          if (face.has_neighbor_ and (cell.global_id_ > face.neighbor_id_))
          {
            if (omega_.Dot(face.normal_) > tolerance)
              std::cout << "  " << cell.global_id_ << " -> "
                        << face.neighbor_id_ << "\n";
            else if (omega_.Dot(face.normal_) < tolerance)
              std::cout << "  " << face.neighbor_id_ << " -> "
                        << cell.global_id_ << "\n";
          } // if outgoing
        }
      }
      std::cout << "\n";
      const auto ghost_ids = grid_.cells.GetGhostGlobalIDs();
      for (const uint64_t global_id : ghost_ids)
      {
        const auto& cell = grid_.cells[global_id];
        for (const auto& face : cell.faces_)
        {
          if (face.has_neighbor_ and (cell.global_id_ > face.neighbor_id_) and
              grid_.IsCellLocal(face.neighbor_id_))
          {
            if (omega_.Dot(face.normal_) > tolerance)
              std::cout << "  " << cell.global_id_ << " -> "
                        << face.neighbor_id_ << "\n";
            else if (omega_.Dot(face.normal_) < tolerance)
              std::cout << "  " << face.neighbor_id_ << " -> "
                        << cell.global_id_ << "\n";
          } // if outgoing
        }
      }
      std::cout << "}\n";

    } // if current location
    Chi::mpi.Barrier();
  } // for p
}