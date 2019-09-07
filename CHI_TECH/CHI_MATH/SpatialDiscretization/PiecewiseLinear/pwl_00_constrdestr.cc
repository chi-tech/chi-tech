#include "pwl.h"
#include "../../Quadratures/quadrature_triangle.h"
#include "../../Quadratures/quadrature_tetrahedron.h"

SpatialDiscretization_PWL::SpatialDiscretization_PWL(int dim)
  : SpatialDiscretization(dim)
{
  CHI_QUADRATURE_TRIANGLE* new_triquad;

  new_triquad = new CHI_QUADRATURE_TRIANGLE(3);
  this->tri_quad_deg5 = new_triquad;

  new_triquad = new CHI_QUADRATURE_TRIANGLE(2,true);
  this->tri_quad_deg3_surf = new_triquad;

  CHI_QUADRATURE_TETRAHEDRON* new_quad;

  new_quad = new CHI_QUADRATURE_TETRAHEDRON(1);
  this->tet_quad_deg1 = new_quad;

  new_quad = new CHI_QUADRATURE_TETRAHEDRON(3);
  this->tet_quad_deg3 = new_quad;

  new_quad = new CHI_QUADRATURE_TETRAHEDRON(3,true);
  this->tet_quad_deg3_surface = new_quad;

  mapping_initialized = false;
}