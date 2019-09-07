#include "quadrature_gausslegendre.h"
#include "LegendrePoly/legendrepoly.h"
#include <math.h>

using namespace chi_math;

void CHI_QUADRATURE_GAUSSLEGENDRE::
    Initialize(int N, int maxiters,
               double tol,bool verbose)
{
  if (verbose)
  {
    printf("Initializing Gauss-Legendre Quadrature with %d q-points\n",N);
  }
  //============================================= InitializeAlphaElements init guess and weights
  double dx = 2.0/2000;
  double xgold = -1.0;
  double fold = Legendre(N,xgold);
  //printf("Fold=%f\n",fold);
  double xgnew = 0.0;
  double fnew = 0.0;
  for (int i=0;i<2000;i++)
  {
    xgnew = xgold + dx;
    fnew = Legendre(N,xgnew);
    //printf("xgnew=%f Fnew=%f  Fold=%f\n",xgnew,fnew,fold);
    if ((fnew*fold)<0.0)
    {
      abscissae.push_back(xgnew);
      weights.push_back(1.0);
    }

    xgold = xgnew;
    fold = fnew;
  }
  //printf("Number of possible roots=%d\n",abscissae.size());

  //============================================= Newton iteration to find root
  for (unsigned k=0; k<abscissae.size(); k++)
  {
    //printf("Finding root %d of %d, ",k+1,N);
    int i=0;
    double xold;
    double xnew;
    double a,b,c;
    double res;
    while (i<maxiters)
    {
      xold = abscissae[k];
      a = Legendre(N,xold);
      b = dLegendredx(N,xold);
      c = 0;
      for (unsigned j=0; j<k; j++)
      {
        c+= 1.0/(xold-abscissae[j]);
      }

      xnew = xold - (a/(b-a*c));

      res = fabs(xnew - xold);
      abscissae[k] = xnew;

      if (res<tol) {break;}
      i++;
    }//while
    weights[k] = 2.0*(1.0-abscissae[k]*abscissae[k])/(N+1)/(N+1)/
                 Legendre(N+1,abscissae[k])/
                 Legendre(N+1,abscissae[k]);

    if (verbose){
      printf("root=%f, weight=%f\n",abscissae[k],weights[k]);
    }

  }
}