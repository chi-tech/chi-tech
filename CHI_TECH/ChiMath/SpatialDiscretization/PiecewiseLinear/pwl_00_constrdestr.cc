#include "pwl.h"
#include "../../Quadratures/quadrature_triangle.h"
#include "../../Quadratures/quadrature_tetrahedron.h"

//###################################################################
/**Constructor.*/
SpatialDiscretization_PWL::
  SpatialDiscretization_PWL(int dim,
                            chi_math::SpatialDiscretizationType sd_method)
  : SpatialDiscretization(dim, sd_method)
{
  chi_math::QuadratureTriangle* new_triquad;

  new_triquad = new chi_math::QuadratureTriangle(3);
  this->tri_quad_deg5 = new_triquad;

  new_triquad = new chi_math::QuadratureTriangle(2,true);
  this->tri_quad_deg3_surf = new_triquad;

  chi_math::QuadratureTetrahedron* new_quad;

  new_quad = new chi_math::QuadratureTetrahedron(1);
  this->tet_quad_deg1 = new_quad;

  new_quad = new chi_math::QuadratureTetrahedron(3);
  this->tet_quad_deg3 = new_quad;

  new_quad = new chi_math::QuadratureTetrahedron(3,true);
  this->tet_quad_deg3_surface = new_quad;

  mapping_initialized = false;
}

