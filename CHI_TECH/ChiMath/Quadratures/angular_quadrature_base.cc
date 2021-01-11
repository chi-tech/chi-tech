#include "angular_quadrature_base.h"

#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include <iomanip>

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

//###################################################################
/**Computes the discrete to moment operator.*/
void chi_math::AngularQuadrature::
  BuildDiscreteToMomentOperator(int scatt_order, bool oneD)
{
  int num_angles = abscissae.size();
  int num_moms = 0;

  d2m_op.clear();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D Slab
  if (oneD)
  {
    int mc=-1; //moment count
    for (int ell=0; ell<=scatt_order; ell++)
    {
      for (int m=0; m<=0; m++)
      {
        std::vector<double> cur_mom; mc++;
        num_moms++;

        for (int n=0; n<num_angles; n++)
        {
          const auto& cur_angle = abscissae[n];
          double value = chi_math::Ylm(ell,m,
                                       cur_angle.phi,
                                       cur_angle.theta);
          double w = weights[n];
          cur_mom.push_back(value*w);
        }

        d2m_op.push_back(cur_mom);
      }//for m
    }//for ell
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D and 3D
  else
  {
    int mc=-1; //moment count
    for (int ell=0; ell<=scatt_order; ell++)
    {
      for (int m=-ell; m<=ell; m++)
      {
        std::vector<double> cur_mom; mc++;
        num_moms++;

        for (int n=0; n<num_angles; n++)
        {
          const auto& cur_angle = abscissae[n];
          double value = chi_math::Ylm(ell,m,
                                       cur_angle.phi,
                                       cur_angle.theta);
          double w = weights[n];
          cur_mom.push_back(value*w);
        }

        d2m_op.push_back(cur_mom);
      }//for m
    }//for ell
  }

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
  BuildMomentToDiscreteOperator(int scatt_order, bool oneD)
{
  int num_angles = abscissae.size();
  int num_moms = 0;

  m2d_op.clear();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D Slab
  if (oneD)
  {
    int mc=-1;
    for (int ell=0; ell<=scatt_order; ell++)
    {
      for (int m=0; m<=0; m++)
      {
        std::vector<double> cur_mom; mc++;
        num_moms++;

        for (int n=0; n<num_angles; n++)
        {
          const auto& cur_angle = abscissae[n];
          double value = ((2.0*ell+1.0)/2.0)*
                         chi_math::Ylm(ell,m,
                                       cur_angle.phi,
                                       cur_angle.theta);
          cur_mom.push_back(value);
        }

        m2d_op.push_back(cur_mom);
      }//for m
    }//for ell
  }
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D and 3D
  else
  {
    int mc=-1;
    for (int ell=0; ell<=scatt_order; ell++)
    {
      for (int m=-ell; m<=ell; m++)
      {
        std::vector<double> cur_mom; mc++;
        num_moms++;

        for (int n=0; n<num_angles; n++)
        {
          const auto& cur_angle = abscissae[n];
          double value = ((2.0*ell+1.0)/4.0/M_PI)*
                         chi_math::Ylm(ell,m,
                                       cur_angle.phi,
                                       cur_angle.theta);
          cur_mom.push_back(value);
        }

        m2d_op.push_back(cur_mom);
      }//for m
    }//for ell
  }

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