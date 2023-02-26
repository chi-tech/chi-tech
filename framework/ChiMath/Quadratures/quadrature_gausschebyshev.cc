#include "quadrature_gausschebyshev.h"
#include <cmath>

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Populates the abscissae and weights for a Gauss-Chebyshev
 * quadrature given the degree \f$ p \f$ of the mononomial such that
 * the quadrature rule integrates exactly the weighted integrand
 * \f$ \rho(x) x^{p} \f$, with \f$ \rho(x) := (1-x^{2})^{-1/2} \f$,
 * on the interval \f$ [-1;+1] \f$.
 * The number of points generated will be ceil((O+1)/2).*/
chi_math::QuadratureGaussChebyshev::
  QuadratureGaussChebyshev(QuadratureOrder in_order,
                           bool verbose)
  : chi_math::Quadrature(in_order)
{
  const unsigned int N = std::ceil(((int)order_ + 1) / 2.0 );
  Initialize(N, verbose);
}

//###################################################################
/**Populates the abscissae and weights for a Gauss-Chebyshev
 * quadrature given the number of desired quadrature points. The
 * order of the quadrature will be 2N-1.*/
chi_math::QuadratureGaussChebyshev::
  QuadratureGaussChebyshev(unsigned int N,
                           bool verbose)
  : chi_math::Quadrature((QuadratureOrder)(2*N-1))
{
  Initialize(N, verbose);
}

//###################################################################
/**Populates the abscissae and weights for a Gauss-Chebyshev
 * quadrature given the number of desired quadrature points.*/
void
chi_math::QuadratureGaussChebyshev::Initialize(unsigned int N, bool verbose)
{
  if (verbose)
    chi::log.Log() << "Initializing Gauss-Chebyshev Quadrature "
                     "with " << N << " q-points";

  const double pi_N = M_PI/N;
  for (unsigned int n = 0; n < N; ++n)
  {
    const double xn = -std::cos( (2*n + 1)*pi_N/2.0 );
    const double wn = pi_N;

    qpoints_.emplace_back(xn);
    weights_.emplace_back(wn);

    if (verbose)
      chi::log.Log()
          << "root[" << n << "]=" << qpoints_[n][0]
          << ", weight=" << weights_[n];
  }

  range_ = {-1, +1};
}
