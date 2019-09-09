#ifndef _diffusion_bndry_vacuum_h
#define _diffusion_bndry_vacuum_h

#include "chi_diffusion_bndry.h"

//###################################################################
/**Robin boundary condition. This type of boundary condition doubles
 * for any boundary condition of the form
 *
\f[
a \phi + b D \hat{n}\cdot \nabla \phi = f
\f]

When \f$ a=0\f$ the boundary condition is equivalent to a <B>Neumann</B>
boundary condition.

\f[
b D\hat{n}\cdot \nabla \phi = f
\f]

 When \f$ a=\frac{1}{4} \f$, \f$ b=\frac{1}{2} \f$
and \f$ f=0 \f$ then the boundary condition is equivalent to a
<B>Vacuum</B> boundary condition.

\f[
\frac{1}{4}\phi + \frac{1}{2}D\hat{n}\cdot \nabla \phi = 0
\f]
 */
class chi_diffusion::BoundaryRobin : public chi_diffusion::Boundary
{
public:
  double a;
  double b;
  double f;
public:
  BoundaryRobin(double a_value, double b_value, double f_value) {
    type = DIFFUSION_ROBIN;
    a = a_value;
    b = b_value;
    f = f_value;
  }
};

#endif