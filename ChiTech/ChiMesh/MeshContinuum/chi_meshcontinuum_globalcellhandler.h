#ifndef CHI_MESHCONTINUUM_GLOBALCELLHANDLER_H_
#define CHI_MESHCONTINUUM_GLOBALCELLHANDLER_H_

#include "ChiMesh/Cell/cell.h"

#include <map>

namespace chi_mesh
{
//##################################################
/**Handles all global index queries.*/
class GlobalCellHandler
{
  friend class MeshContinuum;
private:
  std::vector<uint64_t>& local_cell_glob_indices;

  std::vector<chi_mesh::Cell*>& native_cells;
  std::vector<chi_mesh::Cell*>& foreign_cells;

  std::map<uint64_t,uint64_t>& global_cell_id_to_native_id_map;
  std::map<uint64_t,uint64_t>& global_cell_id_to_foreign_id_map;


private:
  explicit GlobalCellHandler(
    std::vector<uint64_t>& in_local_cell_glob_indices,
    std::vector<chi_mesh::Cell*>& in_native_cells,
    std::vector<chi_mesh::Cell*>& in_foreign_cells,
    std::map<uint64_t,uint64_t>& in_global_cell_id_to_native_id_map,
    std::map<uint64_t,uint64_t>& in_global_cell_id_to_foreign_id_map) :
    local_cell_glob_indices(in_local_cell_glob_indices),
    native_cells(in_native_cells),
    foreign_cells(in_foreign_cells),
    global_cell_id_to_native_id_map(in_global_cell_id_to_native_id_map),
    global_cell_id_to_foreign_id_map(in_global_cell_id_to_foreign_id_map)
  {}

public:
  void push_back(chi_mesh::Cell* new_cell);
  chi_mesh::Cell& operator[](uint64_t cell_global_index);
  const chi_mesh::Cell& operator[](uint64_t cell_global_index) const;

  size_t GetNumGhosts() {return global_cell_id_to_foreign_id_map.size();}

  std::vector<uint64_t> GetGhostGlobalIDs();

  uint64_t GetGhostLocalID(int cell_global_index);
};

}//namespace chi_mesh

#endif //CHI_MESHCONTINUUM_GLOBALCELLHANDLER_H_