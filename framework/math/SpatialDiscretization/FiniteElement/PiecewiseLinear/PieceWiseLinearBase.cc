#include "PieceWiseLinearBase.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/CellMappings/PieceWiseLinear/PieceWiseLinearSlabMapping.h"
#include "math/SpatialDiscretization/CellMappings/PieceWiseLinear/PieceWiseLinearPolygonMapping.h"
#include "math/SpatialDiscretization/CellMappings/PieceWiseLinear/PieceWiseLinearPolyhedronMapping.h"

#include "math/Quadratures/quadrature_line.h"
#include "math/Quadratures/quadrature_triangle.h"
#include "math/Quadratures/quadrature_tetrahedron.h"

#include "math/Quadratures/QuadratureWarehouse.h"

#define UnsupportedCellType(fname)                                             \
  std::invalid_argument((fname) +                                              \
                        ": Unsupported cell type encountered. type_id=" +      \
                        std::to_string(static_cast<int>(cell.Type())));

namespace chi_math::spatial_discretization
{
// ###################################################################
/**Constructor*/
PieceWiseLinearBase::PieceWiseLinearBase(const chi_mesh::MeshContinuum& grid,
                                         QuadratureOrder q_order,
                                         SDMType sdm_type,
                                         CoordinateSystemType cs_type)
  : FiniteElementBase(grid, cs_type, sdm_type, q_order)
{
}

void PieceWiseLinearBase::CreateCellMappings()
{
  constexpr std::string_view fname = __PRETTY_FUNCTION__;

  typedef cell_mapping::PieceWiseLinearSlabMapping SlabSlab;
  typedef cell_mapping::PieceWiseLinearPolygonMapping Polygon;
  typedef cell_mapping::PieceWiseLinearPolyhedronMapping Polyhedron;

  auto& q_warehouse = QuadratureWarehouse::GetInstance();

  auto MakeCellMapping = [this, fname, &q_warehouse](const chi_mesh::Cell& cell)
  {
    using namespace std;
    using namespace chi_math;
    std::unique_ptr<chi_math::CellMapping> mapping;

    switch (cell.Type())
    {
      case chi_mesh::CellType::SLAB:
      {
        const auto& vol_quad =
          q_warehouse.GetQuadrature<chi_math::QuadratureLine>(q_order_,
                                                              {0.0, 1.0});

        mapping = make_unique<SlabSlab>(
          cell, ref_grid_, dynamic_cast<const QuadratureLine&>(vol_quad),
                                coord_sys_type_);
        mapping->Initialize();
        break;
      }
      case chi_mesh::CellType::POLYGON:
      {
        const auto& vol_quad =
          q_warehouse.GetQuadrature<QuadratureTriangle>(q_order_);
        const auto& area_quad =
          q_warehouse.GetQuadrature<chi_math::QuadratureLine>(q_order_,
                                                              {0.0, 1.0});

        mapping = make_unique<Polygon>(
          cell,
          ref_grid_,
          dynamic_cast<const QuadratureTriangle&>(vol_quad),
          dynamic_cast<const QuadratureLine&>(area_quad),
          coord_sys_type_);
        mapping->Initialize();
        break;
      }
      case chi_mesh::CellType::POLYHEDRON:
      {
        const auto& vol_quad =
          q_warehouse.GetQuadrature<QuadratureTetrahedron>(q_order_);
        const auto& area_quad =
          q_warehouse.GetQuadrature<QuadratureTriangle>(q_order_);

        mapping = make_unique<Polyhedron>(
          cell,
          ref_grid_,
          dynamic_cast<const QuadratureTetrahedron&>(vol_quad),
          dynamic_cast<const QuadratureTriangle&>(area_quad),
          coord_sys_type_);
        mapping->Initialize();
        break;
      }
      default:
        throw UnsupportedCellType(std::string(fname))
    }
    return mapping;
  };

  for (const auto& cell : ref_grid_.local_cells)
    cell_mappings_.push_back(MakeCellMapping(cell));

  const auto ghost_ids = ref_grid_.cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto ghost_mapping = MakeCellMapping(ref_grid_.cells[ghost_id]);
    nb_cell_mappings_.insert(
      std::make_pair(ghost_id, std::move(ghost_mapping)));
  }
}
} // namespace chi_math::spatial_discretization
