#include "pwl.h"

//###################################################################
/**Constructor.*/
SpatialDiscretization_PWL::
  SpatialDiscretization_PWL(int dim,
                            chi_math::SpatialDiscretizationType sd_method) :
                            SpatialDiscretization(dim, sd_method),
                            tet_quad_order2(chi_math::QuadratureOrder::SECOND),
                            tet_quad_order2_surface(chi_math::QuadratureOrder::SECOND,true)
{
  chi_math::QuadratureTriangle* new_triquad;

  new_triquad = new chi_math::QuadratureTriangle(3);
  this->tri_quad_deg5 = new_triquad;

  new_triquad = new chi_math::QuadratureTriangle(2,true);
  this->tri_quad_deg3_surf = new_triquad;

//  chi_math::QuadratureTetrahedron* new_quad;

//  new_quad = new chi_math::QuadratureTetrahedron(chi_math::FIRST);
//  this->tet_quad_order1 = new_quad;

//  new_quad = new chi_math::QuadratureTetrahedron(chi_math::SECOND);
//  this->tet_quad_order2 = new_quad;

//  new_quad = new chi_math::QuadratureTetrahedron(chi_math::SECOND,true);
//  this->tet_quad_order2_surface = new_quad;

  mapping_initialized = false;
}

