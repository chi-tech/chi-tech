#include "LagrangeBase.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/CellMappings/Lagrange/LagrangeSlabMapping.h"
#include "math/SpatialDiscretization/CellMappings/Lagrange/LagrangeQuadMapping.h"
#include "math/SpatialDiscretization/CellMappings/Lagrange/LagrangeTriangleMapping.h"
#include "math/SpatialDiscretization/CellMappings/Lagrange/LagrangeHexMapping.h"
#include "math/SpatialDiscretization/CellMappings/Lagrange/LagrangeWedgeMapping.h"
#include "math/SpatialDiscretization/CellMappings/Lagrange/LagrangeTetMapping.h"

#include "math/SpatialDiscretization/CellMappings/PieceWiseLinear/PieceWiseLinearPolygonMapping.h"
#include "math/SpatialDiscretization/CellMappings/PieceWiseLinear/PieceWiseLinearPolyhedronMapping.h"

#include "math/Quadratures/point_quadrature.h"
#include "math/Quadratures/quadrature_line.h"
#include "math/Quadratures/quadrature_triangle.h"
#include "math/Quadratures/quadrature_quadrilateral.h"
#include "math/Quadratures/quadrature_tetrahedron.h"
#include "math/Quadratures/quadrature_hexahedron.h"
#include "math/Quadratures/quadrature_wedge.h"

#include "math/Quadratures/QuadratureWarehouse.h"

#include "chi_log_exceptions.h"

namespace chi_math::spatial_discretization
{

LagrangeBase::LagrangeBase(const chi_mesh::MeshContinuum& grid,
                           QuadratureOrder q_order,
                           SDMType sdm_type,
                           CoordinateSystemType cs_type)
  : FiniteElementBase(grid, cs_type, sdm_type, q_order)
{
  CreateCellMappings();
}

void LagrangeBase::CreateCellMappings()
{
  typedef cell_mapping::LagrangeSlabMapping Slab;
  typedef cell_mapping::LagrangeQuadMapping Quad;
  typedef cell_mapping::LagrangeTriangleMapping Triangle;
  typedef cell_mapping::LagrangeHexMapping Hex;
  typedef cell_mapping::LagrangeWedgeMapping Wedge;
  typedef cell_mapping::LagrangeTetMapping Tetrahedron;

  typedef cell_mapping::PieceWiseLinearPolygonMapping Polygon;
  typedef cell_mapping::PieceWiseLinearPolyhedronMapping Polyhedron;

  auto& q_warehouse = QuadratureWarehouse::GetInstance();

  auto MakeCellMapping = [this, &q_warehouse](const chi_mesh::Cell& cell)
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
                                                              {-1.0, 1.0});
        const auto& area_quad =
          q_warehouse.GetQuadrature<chi_math::PointQuadrature>(q_order_);

        mapping = make_unique<Slab>(
          ref_grid_, cell, vol_quad, area_quad, coord_sys_type_);
        mapping->Initialize();
        break;
      }
      case chi_mesh::CellType::POLYGON:
      {
        if (cell.SubType() == chi_mesh::CellType::QUADRILATERAL)
        {
          const auto& vol_quad =
            q_warehouse.GetQuadrature<QuadratureQuadrilateral>(q_order_);
          const auto& area_quad =
            q_warehouse.GetQuadrature<QuadratureLine>(q_order_, {-1.0, 1.0});

          mapping = make_unique<Quad>(
            ref_grid_, cell, vol_quad, area_quad, coord_sys_type_);
          mapping->Initialize();
          break;
        }
        else if (cell.SubType() == chi_mesh::CellType::TRIANGLE)
        {
          const auto& vol_quad =
            q_warehouse.GetQuadrature<QuadratureTriangle>(q_order_);
          const auto& area_quad =
            q_warehouse.GetQuadrature<QuadratureLine>(q_order_, {-1.0, 1.0});

          mapping = make_unique<Triangle>(
            ref_grid_, cell, vol_quad, area_quad, coord_sys_type_);
          mapping->Initialize();
          break;
        }
        else
        {
          const auto& vol_quad =
            q_warehouse.GetQuadrature<QuadratureTriangle>(q_order_);
          const auto& area_quad =
            q_warehouse.GetQuadrature<QuadratureLine>(q_order_, {-1.0, 1.0});

          mapping = make_unique<Polygon>(
            cell,
            ref_grid_,
            dynamic_cast<const QuadratureTriangle&>(vol_quad),
            dynamic_cast<const QuadratureLine&>(area_quad),
            coord_sys_type_);
          mapping->Initialize();
          break;
        }
      }
      case chi_mesh::CellType::POLYHEDRON:
      {
        if (cell.SubType() == chi_mesh::CellType::HEXAHEDRON)
        {
          const auto& vol_quad =
            q_warehouse.GetQuadrature<QuadratureHexahedron>(q_order_);
          const auto& area_quad =
            q_warehouse.GetQuadrature<QuadratureQuadrilateral>(q_order_);

          mapping = make_unique<Hex>(ref_grid_, cell, vol_quad, area_quad,
                                     coord_sys_type_);
          mapping->Initialize();
          break;
        }
        else if (cell.SubType() == chi_mesh::CellType::WEDGE)
        {
          const auto& vol_quad =
            q_warehouse.GetQuadrature<QuadratureWedge>(q_order_);
          const auto& area_quad1 =
            q_warehouse.GetQuadrature<QuadratureQuadrilateral>(q_order_);
          const auto& area_quad2 =
            q_warehouse.GetQuadrature<QuadratureTriangle>(q_order_);

          mapping = make_unique<Wedge>(
            ref_grid_, cell, vol_quad, area_quad1, area_quad2,
            coord_sys_type_);
          mapping->Initialize();
          break;
        }
        else if (cell.SubType() == chi_mesh::CellType::TETRAHEDRON)
        {
          const auto& vol_quad =
            q_warehouse.GetQuadrature<QuadratureTetrahedron>(q_order_);
          const auto& area_quad =
            q_warehouse.GetQuadrature<QuadratureTriangle>(q_order_);

          mapping =
            make_unique<Tetrahedron>(ref_grid_, cell, vol_quad, area_quad,
                                     coord_sys_type_);
          mapping->Initialize();
          break;
        }
        else
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
      }
      default:
        ChiInvalidArgument("Unsupported cell type encountered");
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