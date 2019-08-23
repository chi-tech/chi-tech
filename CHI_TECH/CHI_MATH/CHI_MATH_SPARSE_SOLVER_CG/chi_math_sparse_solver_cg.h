#ifndef CHI_MATH_SPARSE_SOLVER_CG_H
#define CHI_MATH_SPARSE_SOLVER_CG_H

#include "../chi_math_incdef.h"
#include "../CHI_MATH_MATRIX/chi_math_matrix.h"
#include "../SparseMatrix/chi_math_sparse_matrix_v1.h"

//############################################################################# Class definition
/** Conjugate gradient solver optimized for a sparse matrix.*/
class CHI_MATH_SPARSE_SOLVER_CG
{
public:
	
private:
	double                      tolerance;          ///< Tolerance required for convergence [Default: 1.0e-12]
	bool                        initialZero;        ///< Flag indicating whether initial guess should be zero.
	double                      initialValue;       ///< Set initial Value
	CHI_MATH_MATRIX*			x_kp1;				///< Iteration vectors for next iteration.				
	CHI_MATH_MATRIX*			x_k;                ///< Iteration vector for the current iteration.
	CHI_MATH_MATRIX*			v;
	CHI_MATH_MATRIX*			z;
	CHI_MATH_MATRIX*			Av;
	CHI_MATH_MATRIX*			tv;
	CHI_MATH_MATRIX*			tz;
	CHI_MATH_MATRIX*			rtemp;
	CHI_MATH_MATRIX*			dcv;
	CHI_MATH_MATRIX*			rmult;
	CHI_MATH_MATRIX*			r;
	CHI_MATH_MATRIX*			rmult2;
	CHI_MATH_MATRIX*			r2;
	CHI_MATH_SPARSE_MATRIX*     A;
private:
								
public:
								CHI_MATH_SPARSE_SOLVER_CG();
	void                        InitialGuess(double value);
	void                        SetTolerance(double value);
	void                        Compute(CHI_MATH_SPARSE_MATRIX* Ain);
	int                         Solve(CHI_MATH_MATRIX*& x,CHI_MATH_MATRIX* b);
};





#endif
