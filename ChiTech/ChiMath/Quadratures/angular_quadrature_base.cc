#include "angular_quadrature_base.h"

#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include <iomanip>
#include <numeric>

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

  for (unsigned int i=0; i<Na; i++)
  {
    const auto abscissa =
      chi_math::QuadraturePointPhiTheta(azimuthal[i], polar[i]);

    abscissae.push_back(abscissa);

    const double weight = in_weights[i];
    weights.push_back(weight);
    weight_sum += weight;

    if (verbose)
    {
      char buf[200];
      sprintf(buf,"Varphi=%.2f Theta=%.2f Weight=%.3e\n",
              abscissa.phi*180.0/M_PI,
              abscissa.theta*180.0/M_PI,
              weight);
      ostr << buf;
    }
  }

  //================================================== Create omega list
  for (const auto& qpoint : abscissae)
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

//###################################################################
/**Populates a map of moment m to the Spherical Harmonic indices
 * required.*/
void chi_math::AngularQuadrature::
MakeHarmonicIndices(int scatt_order, int dimension)
{
  if (m_to_ell_em_map.empty())
  {
    if (dimension == 1)
      for (int ell=0; ell<=scatt_order; ell++)
        m_to_ell_em_map.emplace_back(ell,0);
    else if (dimension == 2)
      for (int ell=0; ell<=scatt_order; ell++)
        for (int m=-ell; m<=ell; m+=2)
        {
          if (ell == 0 or m != 0)
            m_to_ell_em_map.emplace_back(ell,m);
        }
    else if (dimension == 3)
      for (int ell=0; ell<=scatt_order; ell++)
        for (int m=-ell; m<=ell; m++)
          m_to_ell_em_map.emplace_back(ell,m);
  }
}

//###################################################################
/**Computes the discrete to moment operator.*/
void chi_math::AngularQuadrature::
  BuildDiscreteToMomentOperator(int scatt_order, int dimension)
{
  if (d2m_op_built) return;

  d2m_op.clear();
  MakeHarmonicIndices(scatt_order,dimension);

  int num_angles = abscissae.size();
  int num_moms = m_to_ell_em_map.size();

  for (const auto& ell_em : m_to_ell_em_map)
  {
    std::vector<double> cur_mom;
    cur_mom.reserve(num_angles);

    for (int n=0; n<num_angles; n++)
    {
      const auto& cur_angle = abscissae[n];
      double value = chi_math::Ylm(ell_em.ell,ell_em.m,
                                   cur_angle.phi,
                                   cur_angle.theta);
      double w = weights[n];
      cur_mom.push_back(value*w);
    }

    d2m_op.push_back(cur_mom);
  }
  d2m_op_built = true;







  //=================================== Verbose printout
  std::stringstream outs;
  outs
    << "\nQuadrature d2m operator:\n";
  for (int n=0; n<num_angles; n++)
  {
    outs << std::setw(5) << n;
    for (int m=0; m<num_moms; m++)
    {
      outs
        << std::setw(15) << std::left << std::fixed
        << std::setprecision(10) << d2m_op[m][n] << " ";
    }
    outs << "\n";
  }

  chi_log.Log(LOG_0VERBOSE_1) << outs.str();
}

//###################################################################
/**Computes the moment to discrete operator.*/
void chi_math::AngularQuadrature::
  BuildMomentToDiscreteOperator(int scatt_order, int dimension)
{
  if (m2d_op_built) return;

  m2d_op.clear();
  MakeHarmonicIndices(scatt_order,dimension);

  int num_angles = abscissae.size();
  int num_moms = m_to_ell_em_map.size();

  const auto normalization =
    std::accumulate(weights.begin(), weights.end(), static_cast<double>(0));

  for (const auto& ell_em : m_to_ell_em_map)
  {
    std::vector<double> cur_mom;
    cur_mom.reserve(num_angles);

    for (int n=0; n<num_angles; n++)
    {
      const auto& cur_angle = abscissae[n];
      double value = ((2.0*ell_em.ell+1.0)/normalization)*
                     chi_math::Ylm(ell_em.ell,ell_em.m,
                                   cur_angle.phi,
                                   cur_angle.theta);
      cur_mom.push_back(value);
    }

    m2d_op.push_back(cur_mom);
  }//for m
  m2d_op_built = true;

  //=================================== Verbose printout
  std::stringstream outs;

  outs
    << "\nQuadrature m2d operator:\n";
  for (int n=0; n<num_angles; n++)
  {
    outs << std::setw(5) << n;
    for (int m=0; m<num_moms; m++)
    {
      outs
        << std::setw(15) << std::left << std::fixed
        << std::setprecision(10) << m2d_op[m][n] << " ";
    }
    outs << "\n";
  }

  chi_log.Log(LOG_0VERBOSE_1) << outs.str();
}