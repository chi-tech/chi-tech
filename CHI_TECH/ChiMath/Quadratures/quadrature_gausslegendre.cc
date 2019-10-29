#include "quadrature_gausslegendre.h"
#include "LegendrePoly/legendrepoly.h"
#include <math.h>

using namespace chi_math;

#include <chi_log.h>

extern ChiLog chi_log;

#include <algorithm>

//###################################################################
/**Initializes the Legendre quadrature.*/
void chi_math::QuadratureGaussLegendre::
    Initialize(int N, int maxiters,
               double tol,bool verbose)
{
  if (verbose)
  {
    printf("Initializing Gauss-Legendre Quadrature with %d q-points\n",N);
  }


  std::vector<double>& roots = abscissae;
  FindRoots(N,roots);

//  std::vector<double> weights2(N);
  weights.resize(N,1.0);
  for (size_t k=0; k<roots.size(); k++)
  {
    weights[k] =
      2.0*(1.0-roots[k]*roots[k])/(N+1)/(N+1)/
      Legendre(N+1,roots[k])/
      Legendre(N+1,roots[k]);

    if (verbose)
      chi_log.Log(LOG_0)
        << "root[" << k << "]=" << roots[k]
        << ", weight=" << weights[k];
  }


//  //============================================= Initialize init guess and weights
//  abscissae.clear();
//  weights.clear();
//  double dx = 2.0/2000;
//  double xgold = -1.0;
//  double fold = Legendre(N,xgold);
//  //printf("Fold=%f\n",fold);
//  double xgnew = 0.0;
//  double fnew = 0.0;
//  for (int i=0;i<2000;i++)
//  {
//    xgnew = xgold + dx;
//    fnew = Legendre(N,xgnew);
//    //printf("xgnew=%f Fnew=%f  Fold=%f\n",xgnew,fnew,fold);
//    if ((fnew*fold)<0.0)
//    {
//      abscissae.push_back(xgnew);
//      weights.push_back(1.0);
//    }
//
//    xgold = xgnew;
//    fold = fnew;
//  }
//  //printf("Number of possible roots=%d\n",abscissae.size());
//
//  //============================================= Newton iteration to find root
//  for (unsigned k=0; k<abscissae.size(); k++)
//  {
//    //printf("Finding root %d of %d, ",k+1,N);
//    int i=0;
//    double xold;
//    double xnew;
//    double a,b,c;
//    double res;
//    while (i<maxiters)
//    {
//      xold = abscissae[k];
//      a = Legendre(N,xold);
//      b = dLegendredx(N,xold);
//      c = 0;
//      for (unsigned j=0; j<k; j++)
//      {
//        c+= 1.0/(xold-abscissae[j]);
//      }
//
//      xnew = xold - (a/(b-a*c));
//
//      res = fabs(xnew - xold);
//      abscissae[k] = xnew;
//
//      if (res<tol) {break;}
//      i++;
//    }//while
//    weights[k] = 2.0*(1.0-abscissae[k]*abscissae[k])/(N+1)/(N+1)/
//                 Legendre(N+1,abscissae[k])/
//                 Legendre(N+1,abscissae[k]);
//
//    if (verbose){
//      printf("root=%f, weight=%f\n",abscissae[k],weights[k]);
//    }
//
//  }
}

//###################################################################
/** Finds the roots of the Legendre polynomial.*/
void chi_math::QuadratureGaussLegendre::FindRoots(
  int N, std::vector<double> &roots)
{
  //======================================== Set overall params
  typedef std::vector<double> Tvecdbl;
  size_t maxiters = 1000;
  double tol = 1.0e-12;
  double adder = 0.99999999*2/std::max(N-1,1);

  //======================================== Populate init guess
  Tvecdbl xn(N, 0.0);
  for(int i=0; i<N; i++)
    xn[i] = -0.99999999 + i*adder;

  //======================================== Find roots
  for (size_t k=0; k<N; k++)
  {
    for (size_t i=0; i<maxiters; i++)
    {
      double xold = xn[k];
      double a = Legendre(N,xold);
      double b = dLegendredx(N,xold);
      double c = 0;

      for (int j=0; j<k; j++)
        c = c+(1.0/(xold-xn[j]));

      double xnew = xold-(a/(b-a*c));
      if (std::isnan(xnew))
      {
        chi_log.Log(LOG_ALLERROR)
          << "QuadratureGaussLegendre::FindRoots "
          << "xnew " << i << " "
          << xnew << " y="
          << a << std::endl;
        exit(EXIT_FAILURE);
      }

      chi_log.Log(LOG_0VERBOSE_2)
        << "xnew " << i << " "
        << xnew << " y="
        << a << std::endl;

      double res = std::fabs(xnew-xold);
      xn[k] = xnew;

      if (res<tol)
        break;

      i = i+1;
    }//while

  }//for k

  std::stable_sort(xn.begin(),xn.end());

  roots = std::move(xn);
}