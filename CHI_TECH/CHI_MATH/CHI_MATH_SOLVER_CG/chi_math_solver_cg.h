#ifndef CHI_MATH_SOLVER_CG_H
#define CHI_MATH_SOLVER_CG_H

#include "../chi_math_incdef.h"
#include "../CHI_MATH_MATRIX/chi_math_matrix.h"

//############################################################################# CLASS DEF
/** Class for efficiently implementing the conjugate gradient solver.
\author CHIV.*/
class CHI_MATH_SOLVER_CG
{
  public:
  		CHI_MATH_MATRIX*  x_kp1;
  		CHI_MATH_MATRIX*  x_k;
  		CHI_MATH_MATRIX*  v;
  		CHI_MATH_MATRIX*  z;
  		CHI_MATH_MATRIX*  Av;
  		CHI_MATH_MATRIX*  tv;
  		CHI_MATH_MATRIX*  tz;
  		CHI_MATH_MATRIX*  rtemp;
  		CHI_MATH_MATRIX*  dcv;
  		CHI_MATH_MATRIX*  rmult;
  		CHI_MATH_MATRIX*  r;
  		CHI_MATH_MATRIX*  rmult2;
  		CHI_MATH_MATRIX*  r2;
  
  public:
	//CHI_MATH_SOLVER_CG(CHI_MATH_MATRIX *x_k);

//00
  ~CHI_MATH_SOLVER_CG();
  //01
  void      Initialize(long int size);
  long int  Solve(CHI_MATH_MATRIX* A, 
                  CHI_MATH_MATRIX*& x, 
                  CHI_MATH_MATRIX* b, 
                  double tolerance=1.0e-12, 
                  bool initialZero=true);
};





#endif
