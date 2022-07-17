#include "ChiMath/Quadratures/cylindrical_angular_quadrature.h"

#include <algorithm>
#include <limits>
#include <numeric>

#include "chi_log.h"

extern ChiLog& chi_log;


chi_math::CylindricalAngularQuadrature::
  CylindricalAngularQuadrature(
    const chi_math::Quadrature& quad_polar,
    const chi_math::Quadrature& quad_azimu,
    const bool verbose)
  : CurvilinearAngularQuadrature()
{
  const auto np = quad_polar.weights.size();
  std::vector<chi_math::Quadrature> quad_azimu_vec(np, quad_azimu);
  Initialize(quad_polar, quad_azimu_vec, verbose);
}


chi_math::CylindricalAngularQuadrature::
  CylindricalAngularQuadrature(
    const chi_math::Quadrature& quad_polar,
    const std::vector<chi_math::Quadrature>& quad_azimu_vec,
    const bool verbose)
  : CurvilinearAngularQuadrature()
{
  Initialize(quad_polar, quad_azimu_vec, verbose);
}


void
chi_math::CylindricalAngularQuadrature::
  Initialize(const chi_math::Quadrature& quad_polar,
             const std::vector<chi_math::Quadrature>& quad_azimu_vec,
             const bool verbose)
{
  //  copies of input quadratures
  auto polar_quad(quad_polar);
  auto azimu_quad_vec(quad_azimu_vec);

  //  --------------------------------------------------------------------------
  //  verifications and corrections (if possible)
  //  --------------------------------------------------------------------------
  const auto eps = std::numeric_limits<double>::epsilon();

  //  consistency among polar quadrature and azimuthal quadratures
  if (polar_quad.weights.size() != azimu_quad_vec.size())
    throw std::invalid_argument("chi_math::CylindricalAngularQuadrature::Initialize : "
                                "number of azimuthal quadratures does not correspond "
                                "to number of polar points of the polar quadrature.");

  //  at present, this class does not handle correctly reduced geometries
  if (polar_quad.weights.size() == 0)
    throw std::invalid_argument("chi_math::CylindricalAngularQuadrature::Initialize : "
                                "invalid polar quadrature size = "
                                +std::to_string(polar_quad.weights.size()));

  for (const auto& azimu_quad : azimu_quad_vec)
    if (azimu_quad.weights.size() == 0)
      throw std::invalid_argument("chi_math::CylindricalAngularQuadrature::Initialize : "
                                  "invalid azimuthal quadrature size = "
                                  +std::to_string(azimu_quad.weights.size()));

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
    throw std::invalid_argument("chi_math::CylindricalAngularQuadrature::Initialize : "
                                "polar quadrature weights sum to zero.");

  //  defined on range [-1;+1]
  if (std::abs(polar_quad.GetRange().first - polar_quad_span.first) > eps ||
      std::abs(polar_quad.GetRange().second- polar_quad_span.second) > eps)
    polar_quad.SetRange(polar_quad_span);

  //  --------------------------------------------------------------------------
  //  verifications on azimuthal quadrature
  //  --------------------------------------------------------------------------
  const double azimu_quad_sum_weights = M_PI;
  const auto azimu_quad_span = std::pair<double, double>(-1, +1);

  for (auto& azimu_quad : azimu_quad_vec)
  {
    //  weights sum to $\pi$
    const auto integral_weights =
      std::accumulate(azimu_quad.weights.begin(), azimu_quad.weights.end(), 0.0);
    if (std::abs(integral_weights) > 0)
    {
      const auto fac = azimu_quad_sum_weights / integral_weights;
      if (std::abs(fac - 1) > eps)
        for (auto& w : azimu_quad.weights)
          w *= fac;
    }
    else
      throw std::invalid_argument("chi_math::CylindricalAngularQuadrature::Initialize : "
                                  "azimuthal quadrature weights sum to zero.");

    //  defined on range [-1;+1]
    if (std::abs(azimu_quad.GetRange().first - azimu_quad_span.first) > eps ||
        std::abs(azimu_quad.GetRange().second- azimu_quad_span.second) > eps)
      azimu_quad.SetRange(azimu_quad_span);

    //  abscissae sorted in ascending order
    auto lt_qp = [](const chi_math::QuadraturePointXYZ& qp0,
                    const chi_math::QuadraturePointXYZ& qp1)
      { return qp0[0] < qp1[0]; };
    if (!std::is_sorted(azimu_quad.qpoints.begin(), azimu_quad.qpoints.end(), lt_qp))
      throw std::invalid_argument("chi_math::CylindricalAngularQuadrature::Initialize : "
                                  "azimuthal quadrature abscissae not in ascending order.");

    //  existence of zero-weight abscissae at the start and at the end of the interval
    if (std::abs(azimu_quad.weights.front()) > eps &&
        std::abs(azimu_quad.qpoints.front()[0] - azimu_quad_span.first) > eps)
    {
      azimu_quad.weights.emplace(azimu_quad.weights.begin(), 0);
      azimu_quad.qpoints.emplace(azimu_quad.qpoints.begin(), azimu_quad_span.first);
    }
    if (std::abs(azimu_quad.weights.back()) > eps &&
        std::abs(azimu_quad.qpoints.back()[0] - azimu_quad_span.second) > eps)
    {
      azimu_quad.weights.emplace(azimu_quad.weights.end(), 0);
      azimu_quad.qpoints.emplace(azimu_quad.qpoints.end(), azimu_quad_span.second);
    }
  }

  //  --------------------------------------------------------------------------
  //  product quadrature : initialisation
  //  --------------------------------------------------------------------------

  //  compute weights, abscissae $(\varphi, \vartheta)$ and direction vectors
  //  $\omega_{pq} := (\mu_{pq}, \xi_{p}, \eta_{pq})$
  weights.clear();
  abscissae.clear();
  omegas.clear();
  for (size_t p = 0; p < azimu_quad_vec.size(); ++p)
  {
    const auto pol_wei = polar_quad.weights[p];
    const auto pol_abs = polar_quad.qpoints[p][0];
    const auto pol_com = std::sqrt(1 - pol_abs * pol_abs);

    for (size_t q = 0; q < azimu_quad_vec[p].weights.size(); ++q)
    {
      const auto& azimu_quad = azimu_quad_vec[p];

      const auto azi_wei = azimu_quad.weights[q];
      const auto azi_abs = azimu_quad.qpoints[q][0];
      const auto azi_com = std::sqrt(1 - azi_abs * azi_abs);

      const auto weight = pol_wei * azi_wei;
      const auto abscissa =
        QuadraturePointPhiTheta(std::acos(azi_abs), std::acos(pol_abs));
      const auto omega =
        chi_mesh::Vector3(pol_com * azi_abs, pol_abs, pol_com * azi_com);

      weights.emplace_back(weight);
      abscissae.emplace_back(abscissa);
      omegas.emplace_back(omega);
    }
  }
  weights.shrink_to_fit();
  abscissae.shrink_to_fit();
  omegas.shrink_to_fit();

  //  map of direction indices
  unsigned int ind0 = 0;
  map_directions.clear();
  for (size_t p = 0; p < azimu_quad_vec.size(); ++p)
  {
    std::vector<unsigned int> vec_directions_p;
    for (size_t q = 0; q < azimu_quad_vec[p].weights.size(); ++q)
      vec_directions_p.emplace_back(ind0 + q);
    map_directions.emplace(p, vec_directions_p);
    ind0 += azimu_quad_vec[p].weights.size();
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
    chi_log.Log(LOG_0) << "curvilinear product quadrature : cylindrical" << std::endl;
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
chi_math::CylindricalAngularQuadrature::InitializeParameters()
{
  const auto pi_sum_q_weights = static_cast<double>(1);

  fac_diamond_difference.resize(weights.size(), 1);
  fac_streaming_operator.resize(weights.size(), 0);
  for (size_t p = 0; p < map_directions.size(); ++p)
  {
    //  interface quantities initialised to starting direction values
    double alpha_interface = 0;
    double phi_interface = abscissae[map_directions[p].front()].phi;
    std::vector<double> mu_interface(2, std::cos(phi_interface));

    //  initialisation permits to forego start direction and final direction
    for (size_t q = 1; q < map_directions[p].size()-1; ++q)
    {
      const auto k = map_directions[p][q];
      const auto w_pq = weights[k];
      const auto mu_pq = omegas[k].x;
      const auto phi_pq = abscissae[k].phi;

      alpha_interface -= w_pq * mu_pq;

      phi_interface -= w_pq * pi_sum_q_weights;
      mu_interface[0] = mu_interface[1];
      mu_interface[1] = std::cos(phi_interface);

      const auto mu = std::cos(phi_pq);
      const auto tau = (mu - mu_interface[0]) / (mu_interface[1] - mu_interface[0]);

      fac_diamond_difference[k] = tau;
      fac_streaming_operator[k] = alpha_interface / (w_pq * tau) + mu_pq;
    }
  }
}


void
chi_math::CylindricalAngularQuadrature::MakeHarmonicIndices(unsigned int scattering_order, int dimension)
{
  if (m_to_ell_em_map.empty())
  {
    if (dimension == 1)
    {
      for (unsigned int l = 0; l <= scattering_order; ++l)
        for (int m = 0; m <= l; ++m)
          if ((l + m) % 2 == 0)
            m_to_ell_em_map.emplace_back(l, m);
    }
    else if (dimension == 2)
      for (unsigned int l = 0; l <= scattering_order; ++l)
        for (int m = 0; m <= l; ++m)
          m_to_ell_em_map.emplace_back(l, m);
    else
      throw std::invalid_argument("chi_math::CylindricalAngularQuadrature::MakeHarmonicIndices : "
                                  "invalid dimension.");
  }
}
