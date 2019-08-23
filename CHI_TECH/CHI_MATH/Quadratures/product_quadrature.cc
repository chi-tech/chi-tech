#include"product_quadrature.h"
#include"quadrature_gausslegendre.h"
#include"quadrature_gausschebyshev.h"

#include <math.h>

//#########################################################
/**Initializes the quadrature with Gauss-Legendre for
 * the polar angles only.*/
void CHI_PRODUCT_QUADRATURE::InitializeWithGL(int Np, bool verbose)
{
  CHI_QUADRATURE_GAUSSLEGENDRE* gl_polar = new CHI_QUADRATURE_GAUSSLEGENDRE;

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
  for (unsigned i=0; i<(1); i++)
  {
    for (unsigned j=0; j<(Np*2); j++)
    {
      QPOINT_PHITHETA* new_pair = new QPOINT_PHITHETA;

      new_pair->phi   = 0.0;
      new_pair->theta = M_PI-acos(gl_polar->abscissae[j]);

      abscissae.push_back(new_pair);

      weight = gl_polar->weights[j];

      weights.push_back(weight);
      weight_sum += weight;

      if (verbose)
      {
        fprintf(stdout,"Varphi=%.2f Theta=%.2f Weight=%.3e\n",
                new_pair->phi*180.0/M_PI,
                new_pair->theta*180.0/M_PI,
                weight);
      }
    }
  }

  //================================================== Create omega list
  for (int n=0; n<abscissae.size(); n++)
  {
    QPOINT_PHITHETA* qpoint = abscissae[n];

    chi_mesh::Vector* new_omega = new chi_mesh::Vector;
    new_omega->x = sin(qpoint->theta)*cos(qpoint->phi);
    new_omega->y = sin(qpoint->theta)*sin(qpoint->phi);
    new_omega->z = cos(qpoint->theta);

    chi_log.Log(LOG_0VERBOSE_1)
    << "Quadrature angle " << n
    << " " << new_omega->PrintS();

    omegas.push_back(new_omega);
  }

  if (verbose) {fprintf(stdout,"Weight sum=%f\n",weight_sum);}

}

//###################################################################
/**Initializes the quadrature with Gauss-Legendre for the polar
 * angles and Gauss-Chebyshev for the azimuthal.*/
void CHI_PRODUCT_QUADRATURE::InitializeWithGC(int Na, bool verbose)
{
  CHI_QUADRATURE_GAUSSCHEBYSHEV* gl_azimu = new CHI_QUADRATURE_GAUSSCHEBYSHEV;

  gl_azimu->Initialize(Na*4);

  double weight     = 0.0;
  double weight_sum = 0.0;

  //================================================= Create azimuthal angles
  for (unsigned i=0; i<(Na*4); i++)
  {
    azimu_ang.push_back(M_PI*(2*(i+1)-1)/(Na*4));
  }

  //================================================== Create polar angles
  polar_ang.push_back(0.5*M_PI);

  //================================================== Create angle pairs
  for (unsigned i=0; i<(Na*4); i++)
  {
    for (unsigned j=0; j<(1); j++)
    {
      QPOINT_PHITHETA* new_pair = new QPOINT_PHITHETA;

      new_pair->phi   = M_PI*(2*(i+1)-1)/(Na*4);
      new_pair->theta = 0.5*M_PI;

      abscissae.push_back(new_pair);

      weight = 2*gl_azimu->weights[i];

      weights.push_back(weight);
      weight_sum += weight;

      if (verbose)
      {
        fprintf(stdout,"Varphi=%.2f Theta=%.2f Weight=%.3e\n",
                new_pair->phi*180.0/M_PI,
                new_pair->theta*180.0/M_PI,
                weight);
      }

    }
  }

  //================================================== Create omega list
  for (int n=0; n<abscissae.size(); n++)
  {
    QPOINT_PHITHETA* qpoint = abscissae[n];

    chi_mesh::Vector* new_omega = new chi_mesh::Vector;
    new_omega->x = sin(qpoint->theta)*cos(qpoint->phi);
    new_omega->y = sin(qpoint->theta)*sin(qpoint->phi);
    new_omega->z = cos(qpoint->theta);

    omegas.push_back(new_omega);
  }

  if (verbose) {fprintf(stdout,"Weight sum=%f\n",weight_sum);}

}

//#########################################################
/**Initializes the quadrature with Gauss-Legendre for
 * both the polar and azimuthal angles.*/
void CHI_PRODUCT_QUADRATURE::InitializeWithGLL(int Na, int Np, bool verbose)
{
  CHI_QUADRATURE_GAUSSLEGENDRE* gl_polar = new CHI_QUADRATURE_GAUSSLEGENDRE;
  CHI_QUADRATURE_GAUSSLEGENDRE* gl_azimu = new CHI_QUADRATURE_GAUSSLEGENDRE;

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
  for (unsigned i=0; i<(Na*4); i++)
  {
    for (unsigned j=0; j<(Np*2); j++)
    {
      QPOINT_PHITHETA* new_pair = new QPOINT_PHITHETA;

      new_pair->phi   = M_PI*gl_azimu->abscissae[i] + M_PI;
      new_pair->theta = M_PI-acos(gl_polar->abscissae[j]);

      abscissae.push_back(new_pair);

      weight = M_PI*gl_azimu->weights[i]*gl_polar->weights[j];

      weights.push_back(weight);
      weight_sum += weight;

      if (verbose)
      {
        fprintf(stdout,"Varphi=%.2f Theta=%.2f Weight=%.3e\n",
                new_pair->phi*180.0/M_PI,
                new_pair->theta*180.0/M_PI,
                weight);
      }
    }
  }

  //================================================== Create omega list
  for (int n=0; n<abscissae.size(); n++)
  {
    QPOINT_PHITHETA* qpoint = abscissae[n];

    chi_mesh::Vector* new_omega = new chi_mesh::Vector;
    new_omega->x = sin(qpoint->theta)*cos(qpoint->phi);
    new_omega->y = sin(qpoint->theta)*sin(qpoint->phi);
    new_omega->z = cos(qpoint->theta);

    omegas.push_back(new_omega);
  }

  if (verbose) {fprintf(stdout,"Weight sum=%f\n",weight_sum);}

}


//###################################################################
/**Initializes the quadrature with Gauss-Legendre for the polar
 * angles and Gauss-Chebyshev for the azimuthal.*/
void CHI_PRODUCT_QUADRATURE::InitializeWithGLC(int Na, int Np, bool verbose)
{
  CHI_QUADRATURE_GAUSSLEGENDRE* gl_polar = new CHI_QUADRATURE_GAUSSLEGENDRE;
  CHI_QUADRATURE_GAUSSCHEBYSHEV* gl_azimu = new CHI_QUADRATURE_GAUSSCHEBYSHEV;

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
  for (unsigned i=0; i<(Na*4); i++)
  {
    for (unsigned j=0; j<(Np*2); j++)
    {
      QPOINT_PHITHETA* new_pair = new QPOINT_PHITHETA;

      new_pair->phi   = M_PI*(2*(i+1)-1)/(Na*4);
      new_pair->theta = M_PI-acos(gl_polar->abscissae[j]);

      abscissae.push_back(new_pair);

      weight = 2*gl_azimu->weights[i]*gl_polar->weights[j];

      weights.push_back(weight);
      weight_sum += weight;

      if (verbose)
      {
        fprintf(stdout,"Varphi=%.2f Theta=%.2f Weight=%.3e\n",
                new_pair->phi*180.0/M_PI,
                new_pair->theta*180.0/M_PI,
                weight);
      }

    }
  }

  //================================================== Create omega list
  for (int n=0; n<abscissae.size(); n++)
  {
    QPOINT_PHITHETA* qpoint = abscissae[n];

    chi_mesh::Vector* new_omega = new chi_mesh::Vector;
    new_omega->x = sin(qpoint->theta)*cos(qpoint->phi);
    new_omega->y = sin(qpoint->theta)*sin(qpoint->phi);
    new_omega->z = cos(qpoint->theta);

    omegas.push_back(new_omega);
  }

  if (verbose) {fprintf(stdout,"Weight sum=%f\n",weight_sum);}

}