#include "fv.h"

#include "CellViews/fv_slab.h"
#include "CellViews/fv_polygon.h"
#include "CellViews/fv_polyhedron.h"

void chi_math::SpatialDiscretization_FV::CreateCellMappings()
{
  constexpr std::string_view fname = "chi_math::SpatialDiscretization_FV::"
                                     "CreateCellMappings";

  auto MakeCellMapping = [this, fname](const chi_mesh::Cell& cell)
  {
    using namespace std;
    using namespace chi_math;
    std::unique_ptr<chi_math::CellMapping> mapping;

    switch (cell.Type())
    {
      case chi_mesh::CellType::SLAB:
        mapping = make_unique<SlabFVValues>(cell, ref_grid);
        break;
      case chi_mesh::CellType::POLYGON:
        mapping = make_unique<PolygonFVValues>(cell, ref_grid);
        break;
      case chi_mesh::CellType::POLYHEDRON:
        mapping = make_unique<PolyhedronFVValues>(cell, ref_grid);
        break;
      default:
        throw std::logic_error(std::string(fname) +
        std::string(": Invalid cell type encountered."));
    }
    return mapping;
  };

  for (const auto& cell : ref_grid->local_cells)
    a_cell_mappings.push_back(MakeCellMapping(cell));

  const auto ghost_ids = ref_grid->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto ghost_mapping = MakeCellMapping(ref_grid->cells[ghost_id]);
    a_nb_cell_mappings.insert(std::make_pair(ghost_id, std::move(ghost_mapping)));
  }
}