#include "quadrature_gausschebyshev.h"
#include <math.h>

void chi_math::QuadratureGaussChebyshev::Initialize(int N,bool verbose)
{
  if (verbose)
  {
    printf("Initializing Gauss-Quadrature Quadrature with %d q-points\n",N);
  }
  for (int n=0;n<N;n++)
  {
    int ns = n + 1;
    double xn = cos( (2.0*ns - 1.0)*M_PI/2.0/N  );
    double wn = M_PI/N;

    if (verbose){
      printf(" root=%f, weight=%f\n",xn,wn);
    }


    abscissae.push_back(xn);
    weights.push_back(wn);
  }
}