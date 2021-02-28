#include "sldfe_sq.h"

#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Performs a simple Riemann integral of a base functor.*/
double chi_math::SimplifiedLDFESQ::Quadrature::RiemannIntegral(BaseFunctor *F, int Ni)
{
  double dangle = M_PI_2/Ni;
  double dtheta = dangle;
  double dphi = dangle;
  double I_riemann = 0.0;
  for (int i=0; i<Ni; ++i)
  {
    double theta = (0.5+i)*dtheta;
    for (int j=0; j<Ni; ++j)
    {
      double phi = (0.5+j)*dphi;
      double mu_r = cos(phi)*sin(theta);
      double eta_r= sin(phi)*sin(theta);
      double xi_r = cos(theta);

      double fval = (*F)(mu_r,eta_r,xi_r);

      I_riemann += fval*sin(theta)*dtheta*dphi;
    }//for j
  }//for i

  return I_riemann;
}

//###################################################################
/**Performs a quadrature integral of a base functor using the
 * supplied SQs.*/
double chi_math::SimplifiedLDFESQ::Quadrature::
QuadratureSSIntegral(BaseFunctor *F)
{
  double I_quadrature = 0.0;
  for (const auto& sq : initial_octant_SQs)
    for (int i=0; i<4; ++i)
    {
      double mu = sq.sub_sqr_points[i][2];
      double theta = acos(mu);
      double phi = acos(sq.sub_sqr_points[i][0]/sin(theta));

      double mu_r = cos(phi)*sin(theta);
      double eta_r= sin(phi)*sin(theta);
      double xi_r = cos(theta);

      double fval = (*F)(mu_r,eta_r,xi_r);

      I_quadrature += sq.sub_sqr_weights[i]*fval;
    }

  return I_quadrature;
}

//###################################################################
/**Performs a test integration of predefined cases.*/
void chi_math::SimplifiedLDFESQ::Quadrature::
TestIntegration(int test_case, double ref_solution, int RiemannN)
{
  struct Case1 : public BaseFunctor
  {double operator()(double mu, double eta, double xi) override
    {return pow(mu,1.0)*pow(eta,1.0)*pow(xi,0.0);}};

  struct Case2 : public BaseFunctor
  {double operator()(double mu, double eta, double xi) override
    {return pow(mu,3.0)*pow(eta,1.0)*pow(xi,1.0);}};

  struct Case3 : public BaseFunctor
  {double operator()(double mu, double eta, double xi) override
    {return pow(mu,3.0)*pow(eta,6.0)*pow(xi,15.0);}};

  struct SphericalHarmonicF : public BaseFunctor
  {double operator()(double mu, double eta, double xi) override
    {
      double theta = acos(xi);
      double phi = acos(mu/sin(theta));

      return chi_math::Ylm(15,3,phi,theta);
    }};

  Case1 case1;
  Case2 case2;
  Case3 case3;
  SphericalHarmonicF SphF;
  BaseFunctor* F = &case1;
  switch (test_case)
  {
    case 1: F = &case1;break;
    case 2: F = &case2;break;
    case 3: F = &case3;break;
    case 4: F = &SphF;break;
    default: F = &case1;
  }

  const int Nd = initial_octant_SQs.size() * 4;
  const int NR = RiemannN;

  double h = 1.0/sqrt(8.0*Nd);
  double I_riemann = ref_solution;
  if (NR>0)
    I_riemann = std::fabs(RiemannIntegral(F,NR));

  double I_quadrature = std::fabs(QuadratureSSIntegral(F));

  char buff0[200],buff1[200],buff2[200];
  sprintf(buff0,"Riemann integral: %.20e\n",I_riemann);
  sprintf(buff1,"Quadrature integral: %.10e\n",I_quadrature);
  sprintf(buff2,"Error_RQ%05d_%06d: %2d %f %e\n",Nd,Nd*8,
         initial_level,h,std::fabs((I_riemann - I_quadrature)/ref_solution));

  chi_log.Log(LOG_0) << buff0;
  chi_log.Log(LOG_0) << buff1;
  chi_log.Log(LOG_0) << buff2;
}