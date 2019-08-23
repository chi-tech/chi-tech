#ifndef _legendrepoly_h
#define _legendrepoly_h

//###################################################################
/**Container object for functions relating to Legendre polynomials.*/
namespace CHI_LEGENDRE_POLYNOMIALS
{
  double Legendre(int N, double x);
  double dLegendredx(int N, double x);

  double AssocLegendre(int ell, int m, double x);

  double Ylm(int ell, int m, double varphi, double theta);
};


#endif