#include "angular_quadrature_base.h"

#include "chi_log.h"

extern ChiLog& chi_log;

//###################################################################
/**Initializes the quadrature with custom angles and weights.*/
void chi_math::AngularQuadrature::
  InitializeWithCustom(std::vector<double> &azimuthal,
                       std::vector<double> &polar,
                       std::vector<double> &in_weights, bool verbose)
{
  size_t Na = azimuthal.size();
  size_t Np = polar.size();
  size_t Nw = in_weights.size();

  if ((Na-Np != 0) or (Na-Nw != 0))
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_math::AngularQuadrature::InitializeWithCustom: supplied"
         " vectors need to be of equal length.";
    exit(EXIT_FAILURE);
  }

  //================================================== Create angle pairs
  std::stringstream ostr;
  double weight_sum = 0.0;

  for (unsigned i=0; i<Na; i++)
  {
    chi_math::QuadraturePointPhiTheta new_pair;

    new_pair.phi   = azimuthal[i];
    new_pair.theta = polar[i];

    abscissae.push_back(new_pair);

    double weight = in_weights[i];
    weights.push_back(weight);
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

  //================================================== Create omega list
  for (auto qpoint : abscissae)
  {
    chi_mesh::Vector3 new_omega;
    new_omega.x = sin(qpoint.theta)*cos(qpoint.phi);
    new_omega.y = sin(qpoint.theta)*sin(qpoint.phi);
    new_omega.z = cos(qpoint.theta);

    omegas.push_back(new_omega);
  }

  if (verbose)
  {
    chi_log.Log(LOG_0)
      << ostr.str() << "\n"
      << "Weight sum=" << weight_sum;
  }

}