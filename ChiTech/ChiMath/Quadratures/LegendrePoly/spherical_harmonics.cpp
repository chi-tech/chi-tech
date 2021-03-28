#include "legendrepoly.h"

#include <cmath>

//###################################################################
/**Implementation of the tesseral spherical harmonics.
 *
 * This code has a whitepaper associated with it
 * <a href="SphericalHarmonics.pdf" target="_blank"><b>Spherical Harmonics</b></a>
 * */
double chi_math::Ylm(unsigned int ell, int m, double varphi, double theta)
{
  auto fac = [](unsigned int x)
  {
    if (x==0) return 1;
    if (x==1) return 1;

    int factorial = 1;
    for (int i=2; i<=x; ++i)
      factorial *= i;

    return factorial;
  };

  if (m<0)
  {
    return sqrt( ( 2.0 ) * fac(ell-std::abs(m))/fac(ell+std::abs(m)) )*
           AssocLegendre(ell,std::abs(m),cos(theta))*
           sin(std::fabs(m)*varphi);
  }
  else if(m==0)
  {
    return AssocLegendre(ell,m,cos(theta));
  }
  else
  {
    return  sqrt( ( 2.0 ) * fac(ell-m)/fac(ell+m) )*
            AssocLegendre(ell,m,cos(theta))*
            cos(m*varphi);
  }

  return 0.0;
}