#ifndef _product_quadrature_h
#define _product_quadrature_h

#include <vector>
#include "../../CHI_MESH/chi_mesh.h"

#define GAUSS_LEGENDRE           1
#define GAUSS_CHEBYSHEV          2
#define GAUSS_LEGENDRE_LEGENDRE  3
#define GAUSS_LEGENDRE_CHEBYSHEV 4

struct QPOINT_PHITHETA
{
  double phi;
  double theta;
};

//######################################################### Class def
/** Parent class for product quadratures*/
class CHI_PRODUCT_QUADRATURE
{
public:
  std::vector<QPOINT_PHITHETA*> abscissae;
  std::vector<double>           polar_ang;
  std::vector<double>           azimu_ang;
  std::vector<double>           weights;
  std::vector<chi_mesh::Vector*>          omegas;

public:
  //product_quadrature.cc
  void InitializeWithGL(int Np, bool verbose=false);
  void InitializeWithGC(int Na, bool verbose=false);
  void InitializeWithGLL(int Na, int Np, bool verbose=false);
  void InitializeWithGLC(int Na, int Np, bool verbose=false);
  int  GetAngleNum(int polar_angle_index, int azimu_ang_index)
  {
    return azimu_ang_index*polar_ang.size() +
           polar_angle_index;
  }

};




#endif