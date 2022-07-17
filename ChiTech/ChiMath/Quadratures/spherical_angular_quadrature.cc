#include "ChiMath/Quadratures/spherical_angular_quadrature.h"

#include <algorithm>
#include <limits>
#include <numeric>

#include "chi_log.h"

extern ChiLog& chi_log;


chi_math::SphericalAngularQuadrature::
  SphericalAngularQuadrature(
    const chi_math::Quadrature& quad_polar,
    const bool verbose)
  : CurvilinearAngularQuadrature()
{
  Initialize(quad_polar, verbose);
}


void
chi_math::SphericalAngularQuadrature::
  Initialize(const chi_math::Quadrature& quad_polar,
             const bool verbose)
{
  //  copies of input quadratures
  auto polar_quad(quad_polar);

  //  --------------------------------------------------------------------------
  //  verifications and corrections (if possible)
  //  --------------------------------------------------------------------------
  const auto eps = std::numeric_limits<double>::epsilon();

  if (polar_quad.weights.size() == 0)
    throw std::invalid_argument("chi_math::SphericalAngularQuadrature::Initialize : "
                                "invalid polar quadrature size = "
                                +std::to_string(polar_quad.weights.size()));

  //  --------------------------------------------------------------------------
  //  verifications on polar quadrature
  //  --------------------------------------------------------------------------
  const double polar_quad_sum_weights = 2;
  const auto polar_quad_span = std::pair<double, double>(-1, +1);

  //  weights sum to 2
  const auto integral_weights =
    std::accumulate(polar_quad.weights.begin(), polar_quad.weights.end(), 0.0);
  if (std::abs(integral_weights) > 0)
  {
    const auto fac = polar_quad_sum_weights / integral_weights;
    if (std::abs(fac - 1) > eps)
      for (auto& w : polar_quad.weights)
        w *= fac;
  }
  else
    throw std::invalid_argument("chi_math::SphericalAngularQuadrature::Initialize : "
                                "polar quadrature weights sum to zero.");

  //  defined on range [-1;+1]
  if (std::abs(polar_quad.GetRange().first - polar_quad_span.first) > eps ||
      std::abs(polar_quad.GetRange().second- polar_quad_span.second) > eps)
    polar_quad.SetRange(polar_quad_span);

  //  abscissae sorted in ascending order
  auto lt_qp = [](const chi_math::QuadraturePointXYZ& qp0,
                  const chi_math::QuadraturePointXYZ& qp1)
    { return qp0[0] < qp1[0]; };
  if (!std::is_sorted(polar_quad.qpoints.begin(), polar_quad.qpoints.end(), lt_qp))
    throw std::invalid_argument("chi_math::SphericalAngularQuadrature::Initialize : "
                                "polar quadrature abscissae not in ascending order.");

  //  existence of zero-weight abscissae at the start and at the end of the interval
  if (std::abs(polar_quad.weights.front()) > eps &&
      std::abs(polar_quad.qpoints.front()[0] - polar_quad_span.first) > eps)
  {
    polar_quad.weights.emplace(polar_quad.weights.begin(), 0);
    polar_quad.qpoints.emplace(polar_quad.qpoints.begin(), polar_quad_span.first);
  }
  if (std::abs(polar_quad.weights.back()) > eps &&
      std::abs(polar_quad.qpoints.back()[0] - polar_quad_span.second) > eps)
  {
    polar_quad.weights.emplace(polar_quad.weights.end(), 0);
    polar_quad.qpoints.emplace(polar_quad.qpoints.end(), polar_quad_span.second);
  }

  //  --------------------------------------------------------------------------
  //  product quadrature : initialisation
  //  --------------------------------------------------------------------------

  //  compute weights, abscissae $(0, \vartheta_{p})$ and direction vectors
  //  $\omega_{p} := ((1-\mu_{p}^{2})^{1/2}, 0, \mu_{p})$
  weights.clear();
  abscissae.clear();
  omegas.clear();
  for (size_t p = 0; p < polar_quad.weights.size(); ++p)
  {
    const auto pol_wei = polar_quad.weights[p];
    const auto pol_abs = polar_quad.qpoints[p][0];
    const auto pol_com = std::sqrt(1 - pol_abs * pol_abs);

    const auto weight = pol_wei;
    const auto abscissa = QuadraturePointPhiTheta(0, std::acos(pol_abs));
    const auto omega = chi_mesh::Vector3(pol_com, 0, pol_abs);

    weights.emplace_back(weight);
    abscissae.emplace_back(abscissa);
    omegas.emplace_back(omega);
  }
  weights.shrink_to_fit();
  abscissae.shrink_to_fit();
  omegas.shrink_to_fit();

  //  map of direction indices
  map_directions.clear();
  for (size_t p = 0; p < polar_quad.weights.size(); ++p)
  {
    std::vector<unsigned int> vec_directions_p;
    vec_directions_p.emplace_back(p);
    map_directions.emplace(p, vec_directions_p);
  }

  //  --------------------------------------------------------------------------
  //  curvilinear product quadrature : compute additional parametrising factors
  //  --------------------------------------------------------------------------
  InitializeParameters();

  //  --------------------------------------------------------------------------
  //  print
  //  --------------------------------------------------------------------------
  if (verbose)
  {
    chi_log.Log(LOG_0) << "map_directions" << std::endl;
    for (const auto& dir : map_directions)
    {
      chi_log.Log(LOG_0) << "polar level " << dir.first << " : ";
      for (const auto& q : dir.second)
        chi_log.Log(LOG_0) << q << ", ";
      chi_log.Log(LOG_0) << std::endl;
    }
    chi_log.Log(LOG_0) << "curvilinear product quadrature : spherical" << std::endl;
    for (size_t k = 0; k < weights.size(); ++k)
      chi_log.Log(LOG_0)
        << "angle index " << k << ": weight = " << weights[k]
        << ", (phi, theta) = (" << abscissae[k].phi << ", " << abscissae[k].theta << ")"
        << ", omega = " << omegas[k].PrintS()
        << ", fac_diamond_difference = " << fac_diamond_difference[k]
        << ", fac_streaming_operator = " << fac_streaming_operator[k] << std::endl;
    const auto sum_weights =
      std::accumulate(weights.begin(), weights.end(), 0.0);
    chi_log.Log(LOG_0) << "sum(weights) = " << sum_weights << std::endl;
  }
}


void
chi_math::SphericalAngularQuadrature::InitializeParameters()
{
  const auto pi_sum_p_weights = M_PI_2;

  fac_diamond_difference.resize(weights.size(), 1);
  fac_streaming_operator.resize(weights.size(), 0);

  //  interface quantities initialised to starting direction values
  double alpha_interface = 0;
  double theta_interface = abscissae[map_directions[0].front()].theta;
  std::vector<double> mu_interface(2, std::cos(theta_interface));

  //  initialisation permits to forego start direction and final direction
  for (size_t p = 1; p < map_directions.size()-1; ++p)
  {
    const auto k = map_directions[p][0];
    const auto w_p = weights[k];
    const auto mu_p = omegas[k].z;
    const auto theta_p = abscissae[k].theta;

    alpha_interface -= w_p * mu_p;

    theta_interface -= w_p * pi_sum_p_weights;
    mu_interface[0] = mu_interface[1];
    mu_interface[1] = std::cos(theta_interface);

    const auto mu = std::cos(theta_p);
    const auto tau = (mu - mu_interface[0]) / (mu_interface[1] - mu_interface[0]);

    fac_diamond_difference[k] = tau;
    fac_streaming_operator[k] = alpha_interface / (w_p * tau) + mu_p;
    fac_streaming_operator[k] *= 2;
  }
}


void
chi_math::SphericalAngularQuadrature::MakeHarmonicIndices(unsigned int scattering_order, int dimension)
{
  if (m_to_ell_em_map.empty())
  {
    if (dimension == 1)
      for (unsigned int l = 0; l <= scattering_order; ++l)
        m_to_ell_em_map.emplace_back(l, 0);
    else
      throw std::invalid_argument("chi_math::SphericalAngularQuadrature::MakeHarmonicIndices : "
                                  "invalid dimension.");
  }
}
