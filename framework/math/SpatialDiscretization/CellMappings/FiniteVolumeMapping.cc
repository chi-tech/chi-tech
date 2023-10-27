#include "math/SpatialDiscretization/CellMappings/FiniteVolumeMapping.h"

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

#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

#include "math/Quadratures/QuadratureWarehouse.h"

#include "chi_log_exceptions.h"

namespace chi_math::cell_mapping
{

FiniteVolumeMapping::FiniteVolumeMapping(
  const chi_mesh::MeshContinuum& grid,
  const chi_mesh::Cell& cell,
  std::vector<std::vector<int>> face_node_mappings,
  CoordinateSystemType coordinate_system_type)
  : CellMapping(grid,
                cell,
                1,
                {cell.centroid_},
                std::move(face_node_mappings),
                coordinate_system_type)
{
}

finite_element::VolumetricQuadraturePointData
FiniteVolumeMapping::MakeVolumetricQuadraturePointData() const
{
  // We divide by the swf purely to be able to reuse logic when
  // computing volume integrals
  auto swf = GetSpatialWeightingFunction();

  return finite_element::VolumetricQuadraturePointData(
    {0},
    {{cell_.centroid_}},
    {{1.0}},
    {{chi_mesh::Vector3(0, 0, 0)}},
    {volume_/swf(cell_.centroid_)},
    face_node_mappings_,
    num_nodes_);
}

finite_element::SurfaceQuadraturePointData
FiniteVolumeMapping::MakeSurfaceQuadraturePointData(size_t face_index) const
{
  // We divide by the swf purely to be able to reuse logic when
  // computing area integrals
  auto swf = GetSpatialWeightingFunction();

  return finite_element::SurfaceQuadraturePointData(
    {0},
    {{chi_mesh::Vector3(0, 0, 0)}},
    {{1.0}},
    {{chi_mesh::Vector3(0, 0, 0)}},
    {areas_[face_index]/swf(cell_.faces_.at(face_index).centroid_)},
    {{chi_mesh::Vector3(0, 0, 0)}},
    face_node_mappings_,
    1);
}

/**This method has been specialized for FV mappings primarily because
* computing the volumes/areas in cylindrical coordinates are fairly difficult
* without using quadrature-rules. This of course could've meant that we
* essentially duplicated the volume computation of all the lagrange elements
* (and PWL elements for polygons and polyhedra). To overcome this we just go
* ahead and temporarily create Lagrange/PWL element mappings for the FV cells
* and make them create the volumes for us.*/
void FiniteVolumeMapping::ComputeCellVolumeAndAreas(
  const chi_mesh::MeshContinuum& grid,
  const chi_mesh::Cell& cell,
  double& volume,
  std::vector<double>& areas)
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

  QuadratureOrder q_order;
  if (coordinate_system_type_ == CoordinateSystemType::CARTESIAN)
    q_order = QuadratureOrder::SECOND;
  else if (coordinate_system_type_ == CoordinateSystemType::CYLINDRICAL)
    q_order = QuadratureOrder::THIRD;
  else if (coordinate_system_type_ == CoordinateSystemType::SPHERICAL)
    q_order = QuadratureOrder::FOURTH;
  else
    ChiInvalidArgument("Unsupported coordinate system encountered");

  using namespace std;

  std::unique_ptr<CellMapping> mapping;

  switch (cell.Type())
  {
    case chi_mesh::CellType::SLAB:
    {
      const auto& vol_quad =
        q_warehouse.GetQuadrature<chi_math::QuadratureLine>(q_order,
                                                            {-1.0, 1.0});
      const auto& area_quad =
        q_warehouse.GetQuadrature<chi_math::PointQuadrature>(q_order);

      mapping = make_unique<Slab>(ref_grid_, cell, vol_quad, area_quad,
                                  coordinate_system_type_);
      mapping->Initialize();
      break;
    }
    case chi_mesh::CellType::POLYGON:
    {
      if (cell.SubType() == chi_mesh::CellType::QUADRILATERAL)
      {
        const auto& vol_quad =
          q_warehouse.GetQuadrature<QuadratureQuadrilateral>(q_order);
        const auto& area_quad =
          q_warehouse.GetQuadrature<QuadratureLine>(q_order, {-1.0, 1.0});

        mapping = make_unique<Quad>(ref_grid_, cell, vol_quad, area_quad,
                                    coordinate_system_type_);
        mapping->Initialize();
        break;
      }
      else if (cell.SubType() == chi_mesh::CellType::TRIANGLE)
      {
        const auto& vol_quad =
          q_warehouse.GetQuadrature<QuadratureTriangle>(q_order);
        const auto& area_quad =
          q_warehouse.GetQuadrature<QuadratureLine>(q_order, {-1.0, 1.0});

        mapping = make_unique<Triangle>(ref_grid_, cell, vol_quad, area_quad,
                                        coordinate_system_type_);
        mapping->Initialize();
        break;
      }
      else
      {
        const auto& vol_quad =
          q_warehouse.GetQuadrature<QuadratureTriangle>(q_order);
        const auto& area_quad =
          q_warehouse.GetQuadrature<QuadratureLine>(q_order, {-1.0, 1.0});

        mapping = make_unique<Polygon>(
          cell,
          ref_grid_,
          dynamic_cast<const QuadratureTriangle&>(vol_quad),
          dynamic_cast<const QuadratureLine&>(area_quad),
          coordinate_system_type_);
        mapping->Initialize();
        break;
      }
    }
    case chi_mesh::CellType::POLYHEDRON:
    {
      if (cell.SubType() == chi_mesh::CellType::HEXAHEDRON)
      {
        const auto& vol_quad =
          q_warehouse.GetQuadrature<QuadratureHexahedron>(q_order);
        const auto& area_quad =
          q_warehouse.GetQuadrature<QuadratureQuadrilateral>(q_order);

        mapping = make_unique<Hex>(ref_grid_, cell, vol_quad, area_quad,
                                   coordinate_system_type_);
        mapping->Initialize();
        break;
      }
      else if (cell.SubType() == chi_mesh::CellType::WEDGE)
      {
        const auto& vol_quad =
          q_warehouse.GetQuadrature<QuadratureWedge>(q_order);
        const auto& area_quad1 =
          q_warehouse.GetQuadrature<QuadratureQuadrilateral>(q_order);
        const auto& area_quad2 =
          q_warehouse.GetQuadrature<QuadratureTriangle>(q_order);

        mapping =
          make_unique<Wedge>(ref_grid_, cell, vol_quad, area_quad1, area_quad2,
                             coordinate_system_type_);
        mapping->Initialize();
        break;
      }
      else if (cell.SubType() == chi_mesh::CellType::TETRAHEDRON)
      {
        const auto& vol_quad =
          q_warehouse.GetQuadrature<QuadratureTetrahedron>(q_order);
        const auto& area_quad =
          q_warehouse.GetQuadrature<QuadratureTriangle>(q_order);

        mapping =
          make_unique<Tetrahedron>(ref_grid_, cell, vol_quad, area_quad,
                                   coordinate_system_type_);
        mapping->Initialize();
        break;
      }
      else
      {
        const auto& vol_quad =
          q_warehouse.GetQuadrature<QuadratureTetrahedron>(q_order);
        const auto& area_quad =
          q_warehouse.GetQuadrature<QuadratureTriangle>(q_order);

        mapping = make_unique<Polyhedron>(
          cell,
          ref_grid_,
          dynamic_cast<const QuadratureTetrahedron&>(vol_quad),
          dynamic_cast<const QuadratureTriangle&>(area_quad),
          coordinate_system_type_);
        mapping->Initialize();
        break;
      }
    }
    default:
      ChiInvalidArgument("Unsupported cell type encountered");
  }

  volume = mapping->CellVolume();
  areas.clear();
  for (size_t f = 0; f < cell.faces_.size(); ++f)
    areas.push_back(mapping->FaceArea(f));
}

} // namespace chi_math::cell_mapping
