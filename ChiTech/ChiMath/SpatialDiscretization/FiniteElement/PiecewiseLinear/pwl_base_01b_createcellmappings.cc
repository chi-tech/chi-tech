#include "pwl_base.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_slab.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_polygon.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_polyhedron.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_slab_cylindrical.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_slab_spherical.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_polygon_cylindrical.h"

#define InvalidCoordinateSystem(fname) \
std::invalid_argument((fname) + \
": Unsupported coordinate system type encountered.");

#define UnsupportedCellType(fname) \
std::invalid_argument((fname) + \
": Unsupported cell type encountered.");

void chi_math::SpatialDiscretization_PWLBase::CreateCellMappings()
{
  constexpr std::string_view fname = "chi_math::SpatialDiscretization_PWLBase::"
                                     "CreateCellMappings";

  typedef SlabMappingFE_PWL                SlabSlab;
  typedef SlabMappingFE_PWL_Cylindrical    SlabCyli;
  typedef SlabMappingFE_PWL_Spherical      SlabSphr;
  typedef PolygonMappingFE_PWL             Polygon;
  typedef PolygonMappingFE_PWL_Cylindrical PolygonCyli;
  typedef PolyhedronMappingFE_PWL          Polyhedron;

  auto MakeCellMapping = [this, fname](const chi_mesh::Cell& cell)
  {
    using namespace std;
    using namespace chi_math;
    std::unique_ptr<chi_math::CellMapping> mapping;

    switch (cell.Type())
    {
      case chi_mesh::CellType::SLAB:
      {
        const auto& vol_quad = line_quad_order_arbitrary;

        switch (cs_type)
        {
          case CoordinateSystemType::CARTESIAN:
            mapping = make_unique<SlabSlab>(cell, ref_grid, vol_quad);
            break;
          case CoordinateSystemType::CYLINDRICAL:
            mapping = make_unique<SlabCyli>(cell, ref_grid, vol_quad);
            break;
          case CoordinateSystemType::SPHERICAL:
            mapping = make_unique<SlabSphr>(cell, ref_grid, vol_quad);
            break;
          default:
            throw InvalidCoordinateSystem(std::string(fname))
        }
        break;
      }
      case chi_mesh::CellType::POLYGON:
      {
        const auto& vol_quad = tri_quad_order_arbitrary;
        const auto& area_quad = line_quad_order_arbitrary;

        switch (cs_type)
        {
          case CoordinateSystemType::CARTESIAN:
            mapping = make_unique<Polygon>(cell, ref_grid, vol_quad, area_quad);
            break;
          case CoordinateSystemType::CYLINDRICAL:
            mapping = make_unique<PolygonCyli>(cell, ref_grid, vol_quad,area_quad);
            break;
          default:
            throw InvalidCoordinateSystem(std::string(fname))
        }
        break;
      }
      case chi_mesh::CellType::POLYHEDRON:
      {
        const auto& vol_quad = tet_quad_order_arbitrary;
        const auto& area_quad = tri_quad_order_arbitrary;

        switch (cs_type)
        {
          case CoordinateSystemType::CARTESIAN:
            mapping = make_unique<Polyhedron>(cell,ref_grid,vol_quad,area_quad);
            break;
          default:
            throw InvalidCoordinateSystem(std::string(fname))
        }
        break;
      }
      default:
        throw UnsupportedCellType(std::string(fname))
    }
    return mapping;
  };

  for (const auto& cell : ref_grid->local_cells)
    cell_mappings.push_back(MakeCellMapping(cell));

  const auto ghost_ids = ref_grid->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto ghost_mapping = MakeCellMapping(ref_grid->cells[ghost_id]);
    nb_cell_mappings.insert(std::make_pair(ghost_id, std::move(ghost_mapping)));
  }
}