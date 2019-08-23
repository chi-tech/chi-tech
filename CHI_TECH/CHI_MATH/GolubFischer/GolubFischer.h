#ifndef _golubfischer_h
#define _golubfischer_h

/*	This module is the encapsulation of the algorithm
 *	depicted in:
 *
 *	[1] Golub G.H. "How to generate unknown orthogonal
 *  	polynomials out of known orthogonal polynomials",
 *  	Numerical Analysis Project, Stanford University,
 *  	November 1991.
 *
 *	Comptuting roots of the polynomial is an adaption of
 *	Newton's method described in:
 *
 *	[2] Barrera-Figueroa V., et al. "Multiple root finder
 *  	algorithm for Legendre and Chebyshev polynomials
 *  	via Newtons method", Annales Mathematicae et
 *  	Informaticae, volume 33, pages 3-13, 2006.
 *
 *	Finally the weights of the resulting Gauss quadrature
 *	is obtained as described in:
 *
 *	[3] Sloan D.P., "A New Multigroup Monte Carlo
 *  	Scattering Algorithm Suitable for Neutral
 *  	and Charged-Particle Boltzmann and
 *  	Fokker-Planck Calculations", SAND83-7094,
 * 		PhD Dissertation, May 1983.
 */

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

typedef std::vector<std::pair<double,double>> Tvecdbl_vecdbl;
typedef std::vector<double> Tvecdbl;



class GolubFischer
{
public:
  Tvecdbl_vecdbl xn_wn;
  Tvecdbl alpha;
  Tvecdbl beta;
public:
  Tvecdbl_vecdbl& GetDiscreteScatAngles(Tvecdbl& mell);
private:
  void MCA(Tvecdbl& in_mell, Tvecdbl& a, Tvecdbl& b, Tvecdbl& c);
  void RootsOrtho(int& N, Tvecdbl& in_alpha, Tvecdbl& in_beta);
  double dOrtho(int ell, double x, Tvecdbl& in_alpha, Tvecdbl& in_beta);
  double Ortho(int ell, double x, Tvecdbl& in_alpha, Tvecdbl& in_beta);
};

#endif