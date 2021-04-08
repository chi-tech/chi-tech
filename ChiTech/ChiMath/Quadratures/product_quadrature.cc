#include "product_quadrature.h"
#include "quadrature_gausslegendre.h"
#include "quadrature_gausschebyshev.h"

#include <cmath>
#include <sstream>

#include <ChiLog/chi_log.h>
extern ChiLog& chi_log;

//#########################################################
/**Initializes the quadrature with Gauss-Legendre for
 * the polar angles only.*/
void chi_math::ProductQuadrature::InitializeWithGL(int Np, bool verbose)
{
  chi_math::QuadratureGaussLegendre gl_polar(Np*2);

  //================================================= Create azimuthal angles
  azimu_ang.clear();
  azimu_ang.emplace_back(0.0);

  //================================================== Create polar angles
  polar_ang.clear();
  for (unsigned int j = 0; j < (Np*2); ++j)
    polar_ang.emplace_back(M_PI-acos(gl_polar.qpoints[j][0]));

  //================================================== Create combined weights
  auto& weights = gl_polar.weights;

  //================================================== Initialize
  InitializeWithCustom(azimu_ang, polar_ang, weights, verbose);
}

//#########################################################
/**Initializes the quadrature with Gauss-Legendre for
 * both the polar and azimuthal angles.*/
void chi_math::ProductQuadrature::InitializeWithGLL(int Na, int Np, bool verbose)
{
  chi_math::QuadratureGaussLegendre gl_polar(Np*2);
  chi_math::QuadratureGaussLegendre gl_azimu(Na*4);

  //================================================= Create azimuthal angles
  azimu_ang.clear();
  for (unsigned int i = 0; i < (Na*4); ++i)
    azimu_ang.emplace_back(M_PI*gl_azimu.qpoints[i][0] + M_PI);

  //================================================== Create polar angles
  polar_ang.clear();
  for (unsigned int j = 0; j < (Np*2); ++j)
    polar_ang.emplace_back(M_PI-acos(gl_polar.qpoints[j][0]));

  //================================================== Create combined weights
  std::vector<double> weights;
  for (unsigned int i = 0; i < azimu_ang.size(); ++i)
    for (unsigned int j = 0; j < polar_ang.size(); ++j)
      weights.emplace_back(M_PI*gl_azimu.weights[i]*gl_polar.weights[j]);

  //================================================== Initialize
  InitializeWithCustom(azimu_ang, polar_ang, weights, verbose);
}

//###################################################################
/**Initializes the quadrature with Gauss-Legendre for the polar
 * angles and Gauss-Chebyshev for the azimuthal.*/
void chi_math::ProductQuadrature::InitializeWithGLC(int Na, int Np, bool verbose)
{
  chi_math::QuadratureGaussLegendre gl_polar(Np*2);
  chi_math::QuadratureGaussChebyshev gc_azimu(Na*4);

  //================================================= Create azimuthal angles
  azimu_ang.clear();
  for (unsigned int i = 0; i < (Na*4); ++i)
    azimu_ang.emplace_back(M_PI*(2*(i+1)-1)/(Na*4));

  //================================================== Create polar angles
  polar_ang.clear();
  for (unsigned int j = 0; j < (Np*2); ++j)
    polar_ang.emplace_back(M_PI-acos(gl_polar.qpoints[j][0]));

  //================================================== Create combined weights
  std::vector<double> weights;
  for (unsigned int i = 0; i < azimu_ang.size(); ++i)
    for (unsigned int j = 0; j < polar_ang.size(); ++j)
      weights.emplace_back(2*gc_azimu.weights[i]*gl_polar.weights[j]);

  //================================================== Initialize
  InitializeWithCustom(azimu_ang, polar_ang, weights, verbose);
}

//###################################################################
/**Initializes the quadrature with custom angles and weights.*/
void chi_math::ProductQuadrature::
  InitializeWithCustom(std::vector<double>& azimuthal,
                       std::vector<double>& polar,
                       std::vector<double>& in_weights, bool verbose)
{
  size_t Na = azimuthal.size();
  size_t Np = polar.size();
  size_t Nw = in_weights.size();

  if (Nw != Na*Np)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Product Quadrature, InitializeWithCustom: mismatch in the amount "
         "angles and weights. Number of azimuthal angles times number "
         "polar angles must equal the amount of weights.";
    exit(EXIT_FAILURE);
  }

  azimu_ang = azimuthal;
  polar_ang = polar;

  if (verbose)
  {
    chi_log.Log(LOG_0) << "Azimuthal angles:";
    for (const auto& ang : azimu_ang)
      chi_log.Log(LOG_0) << ang;

    chi_log.Log(LOG_0) << "Polar angles:";
    for (const auto& ang : polar_ang)
      chi_log.Log(LOG_0) << ang;
  }

  //================================================== Create angle pairs
  map_directions.clear();
  for (unsigned int j = 0; j < Np; ++j)
    map_directions.emplace(j, std::vector<unsigned int>());

  abscissae.clear();
  weights.clear();
  std::stringstream ostr;
  double weight_sum = 0.0;
  for (unsigned int i = 0; i < Na; ++i)
  {
    for (unsigned int j = 0; j < Np; ++j)
    {
      map_directions[j].emplace_back(i*Np+j);

      chi_math::QuadraturePointPhiTheta new_pair;

      new_pair.phi   = azimu_ang[i];
      new_pair.theta = polar_ang[j];

      abscissae.emplace_back(new_pair);

      double weight = in_weights[i*Np+j];
      weights.emplace_back(weight);
      weight_sum += weight;

      if (verbose)
      {
        char buf[200];
        sprintf(buf,"Varphi=%.2f Theta=%.2f Weight=%.3e\n",
                new_pair.phi*180.0/M_PI,
                new_pair.theta*180.0/M_PI,
                weight);
        ostr << buf;
      }
    }
  }

  //================================================== Create omega list
  omegas.clear();
  for (const auto& qpoint : abscissae)
  {
    chi_mesh::Vector3 new_omega;
    new_omega.x = sin(qpoint.theta)*cos(qpoint.phi);
    new_omega.y = sin(qpoint.theta)*sin(qpoint.phi);
    new_omega.z = cos(qpoint.theta);

    omegas.emplace_back(new_omega);

    if (verbose)
      chi_log.Log(LOG_0) << "Quadrature angle=" << new_omega.PrintS();
  }

  if (verbose)
  {
    chi_log.Log(LOG_0)
      << ostr.str() << "\n"
      << "Weight sum=" << weight_sum;
  }

}
