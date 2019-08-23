#ifndef _quadrature_gauss_chebyshev_h
#define _quadrature_gauss_chebyshev_h

#include "quadrature.h"

//######################################################### Class Def
/**Gauss-Chebyshev quadrature.*/
class CHI_QUADRATURE_GAUSSCHEBYSHEV : public CHI_QUADRATURE
{
public:
  //01
  void Initialize(int N, bool verbose=false);

};

#endif