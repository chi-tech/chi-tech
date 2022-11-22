#ifndef SPATIAL_DISCRETIZATION_PWL_BASE_H
#define SPATIAL_DISCRETIZATION_PWL_BASE_H

#include "ChiMath/SpatialDiscretization/FiniteElement/spatial_discretization_FE.h"

#include "ChiMath/Quadratures/quadrature_line.h"
#include "ChiMath/Quadratures/quadrature_triangle.h"
#include "ChiMath/Quadratures/quadrature_quadrilateral.h"
#include "ChiMath/Quadratures/quadrature_tetrahedron.h"

namespace chi_math
{

  class SpatialDiscretization_PWLBase : public chi_math::SpatialDiscretization_FE
  {
  protected:
    QuadratureLine          line_quad_order_arbitrary;
    QuadratureTriangle      tri_quad_order_arbitrary;
    QuadratureQuadrilateral quad_quad_order_arbitrary;
    QuadratureTetrahedron   tet_quad_order_arbitrary;

  protected:
    explicit
    SpatialDiscretization_PWLBase(chi_mesh::MeshContinuumPtr& in_grid,
                                  finite_element::SetupFlags in_setup_flags,
                                  QuadratureOrder qorder,
                                  SDMType in_type,
                                  CoordinateSystemType in_cs_type
                                  );
    //01
  protected:
    void PreComputeCellSDValues();
    void PreComputeNeighborCellSDValues();

    void CreateCellMappings();
  };

}//namespace chi_math

#endif //SPATIAL_DISCRETIZATION_PWL_BASE_H
