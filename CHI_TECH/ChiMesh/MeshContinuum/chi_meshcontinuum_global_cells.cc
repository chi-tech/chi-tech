#include "chi_meshcontinuum.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

//###################################################################
/**Adds a new cell to grid registry.*/
void chi_mesh::MeshContinuum::GlobalCellHandler::
push_back(chi_mesh::Cell *new_cell)
{
  if (new_cell->partition_id == chi_mpi.location_id)
  {
    local_cells.local_cell_ind.push_back(new_cell->global_id);
    int local_cell_index = local_cells.local_cell_ind.size()-1;
    new_cell->local_id = local_cell_index;

    local_cells.native_cells.push_back(new_cell);

    std::pair<int,int> new_info = std::make_pair(
      new_cell->global_id,
      local_cells.native_cells.size()-1);

    global_cell_native_index_set.insert(new_info);
  }
  else
  {
    local_cells.foreign_cells.push_back(new_cell);

    std::pair<int,int> new_info = std::make_pair(
      new_cell->global_id,
      local_cells.foreign_cells.size() - 1);

    global_cell_foreign_index_set.insert(new_info);
  }

}

//###################################################################
/**Returns a pointer to a cell given its global cell index.*/
chi_mesh::Cell* &chi_mesh::MeshContinuum::GlobalCellHandler::
operator[](int cell_global_index)
{
  //======================================== Define functor
  struct Functor
  {
    int cell_g_id;

    Functor(int in_cgi) : cell_g_id(in_cgi) {}

    bool operator()(const std::pair<int,int>& eval)
    { return eval.first == cell_g_id; }
  }functor(cell_global_index);

  //======================================== First look in native cells
  auto low_bound = global_cell_native_index_set.upper_bound(
    std::make_pair(cell_global_index-1,0)
  );

  auto upp_bound = global_cell_native_index_set.upper_bound(
    std::make_pair(cell_global_index+1,0)
  );

  auto native_loc = std::find_if(low_bound,
                                 upp_bound,
                                 functor);

  //======================================== Then look in foreign cells
  if (native_loc == upp_bound)
  {
    auto low_bound = global_cell_foreign_index_set.upper_bound(
      std::make_pair(cell_global_index-1,0)
    );

    auto upp_bound = global_cell_foreign_index_set.upper_bound(
      std::make_pair(cell_global_index+1,0)
    );

    auto loc = std::find_if(low_bound,
                            upp_bound,
                            functor);

    if (loc != upp_bound)
      return local_cells.foreign_cells[loc->second];
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "chi_mesh::MeshContinuum::cells. Mapping error."
        << "\n"
        << cell_global_index;

      exit(EXIT_FAILURE);
    }
  }
  else
    return local_cells.native_cells[native_loc->second];
}

//###################################################################
/**Returns the total number of global cells.*/
size_t chi_mesh::MeshContinuum::GetGlobalNumberOfCells()
{
  int num_local_cells = local_cells.size();
  int num_globl_cells = 0;

  MPI_Allreduce(&num_local_cells,
                &num_globl_cells,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD);

  return num_globl_cells;
}

//###################################################################
/**Get the number of ghost cells. These are cells that
 * neighbors to this partition's cells but are on a different
 * partition.*/
int chi_mesh::MeshContinuum::GlobalCellHandler::
  GetNumGhosts()
{
  return global_cell_foreign_index_set.size();
}

//###################################################################
/**Returns the cell global ids of all ghost cells. These are cells that
 * neighbors to this partition's cells but are on a different
 * partition.*/
std::vector<int> chi_mesh::MeshContinuum::GlobalCellHandler::
  GetGhostGlobalIDs()
{
  std::vector<int> ids;
  ids.reserve(GetNumGhosts());

//  for (auto map : global_cell_foreign_index_set)
//    ids.push_back(map.first);
  for (auto cell : local_cells.foreign_cells)
    ids.push_back(cell->global_id);

  return ids;
}

//###################################################################
/**Returns the local storage address of a ghost cell. If the
 * ghost is not truly a ghost then -1 is returned, but is wasteful and
 * therefore the user of this function should implement code
 * to prevent it.*/
int chi_mesh::MeshContinuum::GlobalCellHandler::
  GetGhostLocalID(int cell_global_index)
{
  int local_id = -1;

  //======================================== Define functor
  struct Functor
  {
    int cell_g_id;

    Functor(int in_cgi) : cell_g_id(in_cgi) {}

    bool operator()(const std::pair<int,int>& eval)
    { return eval.first == cell_g_id; }
  }functor(cell_global_index);

  //======================================== Search foreign cells
  auto low_bound = global_cell_foreign_index_set.upper_bound(
    std::make_pair(cell_global_index-1,0)
  );

  auto upp_bound = global_cell_foreign_index_set.upper_bound(
    std::make_pair(cell_global_index+1,0)
  );

  auto loc = std::find_if(low_bound,
                          upp_bound,
                          functor);

  if (loc != upp_bound)
    local_id = loc->second;
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "Grid GetGhostLocalID failed to find cell " << cell_global_index;
    exit(EXIT_FAILURE);
  }

  return local_id;
}