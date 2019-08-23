#include "legendrepoly.h"
#include <math.h>

#define PI 3.14159265359

double fac(int x)
{
  if (x==0) return 1;
  if (x==1) return 1;

  return fac(x-1)*x;
}

//###################################################################
/**Implementation of the tesseral spherical harmonics.
 *
 * This code has a whitepaper associated with it
 * <a href="SphericalHarmonics.pdf" target="_blank"><b>Spherical Harmonics</b></a>
 * */
double CHI_LEGENDRE_POLYNOMIALS::Ylm(int ell, int m, double varphi, double theta)
{
  double el = 1.0*ell;
  double em = 1.0*m;

  double C = (2*el + 1.0)/4.0/M_PI;
  C=1.0;
  if (em<0)
  {
    return
               sqrt( ( 2.0 ) * fac(el-abs(em))/fac(el+abs(em)) )*
               AssocLegendre(el,abs(em),cos(theta))*
               sin(fabs(em)*varphi);
  }
  else if(em==0)
  {
    return
               sqrt( ( 1.0 ) * fac(el-em)/fac(el+em) )*
               AssocLegendre(el,em,cos(theta));
  }
  else
  {
    return
               sqrt( ( 2.0 ) * fac(el-em)/fac(el+em) )*
               AssocLegendre(el,em,cos(theta))*
               cos(em*varphi);
  }

}