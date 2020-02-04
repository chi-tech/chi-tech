#include"product_quadrature.h"
#include"quadrature_gausslegendre.h"
#include"quadrature_gausschebyshev.h"

#include <cmath>
#include <sstream>

#include <ChiLog/chi_log.h>
extern ChiLog chi_log;

//#########################################################
/**Initializes the quadrature with Gauss-Legendre for
 * the polar angles only.*/
void chi_math::ProductQuadrature::InitializeWithGL(int Np, bool verbose)
{
  auto gl_polar = new chi_math::QuadratureGaussLegendre;

  gl_polar->Initialize(Np*2);

  double weight     = 0.0;
  double weight_sum = 0.0;

  //================================================= Create azimuthal angles
  azimu_ang.push_back(0.0);

  //================================================== Create polar angles
  for (unsigned j=0; j<(Np*2); j++)
  {
    polar_ang.push_back(M_PI-acos(gl_polar->abscissae[j]));
  }

  //================================================== Create angle pairs
  std::stringstream ostr;
  for (unsigned i=0; i<(1); i++)
  {
    for (unsigned j=0; j<(Np*2); j++)
    {
      auto new_pair = new chi_math::QuadraturePointPhiTheta;

      new_pair->phi   = 0.0;
      new_pair->theta = M_PI-acos(gl_polar->abscissae[j]);

      abscissae.push_back(new_pair);

      weight = gl_polar->weights[j];

      weights.push_back(weight);
      weight_sum += weight;

      if (verbose)
      {
        char buf[200];
        sprintf(buf,"Varphi=%.2f Theta=%.2f Weight=%.3e\n",
                new_pair->phi*180.0/M_PI,
                new_pair->theta*180.0/M_PI,
                weight);
        ostr << buf;
      }
    }
  }

  //================================================== Create omega list
  for (size_t n=0; n<abscissae.size(); n++)
  {
    chi_math::QuadraturePointPhiTheta* qpoint = abscissae[n];

    auto new_omega = new chi_mesh::Vector;
    new_omega->x = sin(qpoint->theta)*cos(qpoint->phi);
    new_omega->y = sin(qpoint->theta)*sin(qpoint->phi);
    new_omega->z = cos(qpoint->theta);

    chi_log.Log(LOG_0VERBOSE_1)
    << "Quadrature angle " << n
    << " " << new_omega->PrintS();

    omegas.push_back(new_omega);
  }

  if (verbose)
  {
    chi_log.Log(LOG_0)
      << ostr.str() << "\n"
      << "Weight sum=" << weight_sum;
  }

}

//#########################################################
/**Initializes the quadrature with Gauss-Legendre for
 * both the polar and azimuthal angles.*/
void chi_math::ProductQuadrature::InitializeWithGLL(int Na, int Np, bool verbose)
{
  auto gl_polar = new chi_math::QuadratureGaussLegendre;
  auto gl_azimu = new chi_math::QuadratureGaussLegendre;

  gl_polar->Initialize(Np*2);
  gl_azimu->Initialize(Na*4);

  double weight     = 0.0;
  double weight_sum = 0.0;

  //================================================= Create azimuthal angles
  for (unsigned i=0; i<(Na*4); i++)
  {
    azimu_ang.push_back(M_PI*gl_azimu->abscissae[i] + M_PI);
  }

  //================================================== Create polar angles
  for (unsigned j=0; j<(Np*2); j++)
  {
    polar_ang.push_back(M_PI-acos(gl_polar->abscissae[j]));
  }

  //================================================== Create angle pairs
  std::stringstream ostr;
  for (unsigned i=0; i<(Na*4); i++)
  {
    for (unsigned j=0; j<(Np*2); j++)
    {
      auto new_pair = new chi_math::QuadraturePointPhiTheta;

      new_pair->phi   = M_PI*gl_azimu->abscissae[i] + M_PI;
      new_pair->theta = M_PI-acos(gl_polar->abscissae[j]);

      abscissae.push_back(new_pair);

      weight = M_PI*gl_azimu->weights[i]*gl_polar->weights[j];

      weights.push_back(weight);
      weight_sum += weight;

      if (verbose)
      {
        char buf[200];
        sprintf(buf,"Varphi=%.2f Theta=%.2f Weight=%.3e\n",
                new_pair->phi*180.0/M_PI,
                new_pair->theta*180.0/M_PI,
                weight);
        ostr << buf;
      }
    }
  }

  //================================================== Create omega list
  for (size_t n=0; n<abscissae.size(); n++)
  {
    chi_math::QuadraturePointPhiTheta* qpoint = abscissae[n];

    auto new_omega = new chi_mesh::Vector;
    new_omega->x = sin(qpoint->theta)*cos(qpoint->phi);
    new_omega->y = sin(qpoint->theta)*sin(qpoint->phi);
    new_omega->z = cos(qpoint->theta);

    omegas.push_back(new_omega);
  }

  if (verbose)
  {
    chi_log.Log(LOG_0)
      << ostr.str() << "\n"
      << "Weight sum=" << weight_sum;
  }

}


//###################################################################
/**Initializes the quadrature with Gauss-Legendre for the polar
 * angles and Gauss-Chebyshev for the azimuthal.*/
void chi_math::ProductQuadrature::InitializeWithGLC(int Na, int Np, bool verbose)
{
  auto gl_polar = new chi_math::QuadratureGaussLegendre;
  auto gl_azimu = new chi_math::QuadratureGaussChebyshev;

  gl_polar->Initialize(Np*2);
  gl_azimu->Initialize(Na*4);

  double weight     = 0.0;
  double weight_sum = 0.0;

  //================================================= Create azimuthal angles
  for (unsigned i=0; i<(Na*4); i++)
  {
    azimu_ang.push_back(M_PI*(2*(i+1)-1)/(Na*4));
  }

  //================================================== Create polar angles
  for (unsigned j=0; j<(Np*2); j++)
  {
    polar_ang.push_back(M_PI-acos(gl_polar->abscissae[j]));
  }

  //================================================== Create angle pairs
  std::stringstream ostr;
  for (unsigned i=0; i<(Na*4); i++)
  {
    for (unsigned j=0; j<(Np*2); j++)
    {
      auto new_pair = new chi_math::QuadraturePointPhiTheta;

      new_pair->phi   = M_PI*(2*(i+1)-1)/(Na*4);
      new_pair->theta = M_PI-acos(gl_polar->abscissae[j]);

      abscissae.push_back(new_pair);

      weight = 2*gl_azimu->weights[i]*gl_polar->weights[j];

      weights.push_back(weight);
      weight_sum += weight;

      if (verbose)
      {
        char buf[200];
        sprintf(buf,"Varphi=%.2f Theta=%.2f Weight=%.3e\n",
                new_pair->phi*180.0/M_PI,
                new_pair->theta*180.0/M_PI,
                weight);
        ostr << buf;
      }
    }
  }

  //================================================== Create omega list
  for (size_t n=0; n<abscissae.size(); n++)
  {
    chi_math::QuadraturePointPhiTheta* qpoint = abscissae[n];

    chi_mesh::Vector* new_omega = new chi_mesh::Vector;
    new_omega->x = sin(qpoint->theta)*cos(qpoint->phi);
    new_omega->y = sin(qpoint->theta)*sin(qpoint->phi);
    new_omega->z = cos(qpoint->theta);

    omegas.push_back(new_omega);
  }

  if (verbose)
  {
    chi_log.Log(LOG_0)
      << ostr.str() << "\n"
      << "Weight sum=" << weight_sum;
  }

}

//###################################################################
/**Initializes the quadrature with Gauss-Legendre for the polar
 * angles and Gauss-Chebyshev for the azimuthal.*/
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

  //================================================== Create angle pairs
  std::stringstream ostr;
  double weight_sum = 0.0;
  int nw = -1;
  for (unsigned i=0; i<Na; i++)
  {
    for (unsigned j=0; j<Np; j++)
    {
      ++nw;
      auto new_pair = new chi_math::QuadraturePointPhiTheta;

      new_pair->phi   = azimuthal[i];
      new_pair->theta = polar[j];

      abscissae.push_back(new_pair);

      double weight = in_weights[nw];
      weights.push_back(weight);
      weight_sum += weight;

      if (verbose)
      {
        char buf[200];
        sprintf(buf,"Varphi=%.2f Theta=%.2f Weight=%.3e\n",
                new_pair->phi*180.0/M_PI,
                new_pair->theta*180.0/M_PI,
                weight);
        ostr << buf;
      }
    }
  }

  //================================================== Create omega list
  for (auto qpoint : abscissae)
  {
    auto new_omega = new chi_mesh::Vector;
    new_omega->x = sin(qpoint->theta)*cos(qpoint->phi);
    new_omega->y = sin(qpoint->theta)*sin(qpoint->phi);
    new_omega->z = cos(qpoint->theta);

    omegas.push_back(new_omega);
  }

  if (verbose)
  {
    chi_log.Log(LOG_0)
      << ostr.str() << "\n"
      << "Weight sum=" << weight_sum;
  }

}