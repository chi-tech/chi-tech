#include "chi_math.h"
#include <assert.h>

//############################################################################# GAUSS ELIMINATION
/** Gauss Elimination without pivoting.*/
void CHI_MATH::GaussElimination(std::vector<std::vector<double> > &A,
	                              std::vector<double> &b, int n)
{
	// Forward elimination
	for(int i = 0;i < n-1;++i)
	{
		const std::vector<double>& ai = A[i];
		double bi = b[i];
		double factor = 1.0/A[i][i];
		for(int j = i+1;j < n;++j)
		{
			std::vector<double>& aj = A[j];
			double val = aj[i] * factor;
			b[j] -= val * bi;
			for(int k = i+1;k < n;++k)
				aj[k] -= val * ai[k];
		}
	}

	// Back substitution
	for(int i = n-1;i >= 0;--i)
	{
		const std::vector<double>& ai = A[i];
		double bi = b[i];
		for(int j = i+1;j < n;++j)
			bi -= ai[j] * b[j];
		b[i] = bi/ai[i];
	}
}

//############################################################################# TDMA
/** Tri-Diagonal Matrix Algorithm solver by direct inversion.
\author Jan*/
void CHI_MATH::TDMA(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* x, CHI_MATH_MATRIX* b, int N)
{
	A->sij(1,3,-A->ij(1,3)/A->ij(1,2));
	b->sij(1,1,b->ij(1,1)/A->ij(1,2));

	for (int i=2;i<=N;i++)
	{
		A->sij(i,3,-A->ij(i,3)/( A->ij(i,2) + A->ij(i,1)*A->ij(i-1,3)    ));
		b->sij(i,1, (b->ij(i,1)-A->ij(i,1)*b->ij(i-1,1))/ ( A->ij(i,2) + A->ij(i,1)*A->ij(i-1,3)    )  );
	}

	for (int i=(N-1);i>=1;i--)
	{
		b->sij(i,1,A->ij(i,3)*b->ij(i+1)+b->ij(i,1));
	}

	for (int i=1;i<=N;i++)
	{
		x->sij(i,1,b->ij(i,1));
	}
}

//############################################################################# Jacobi solver
/** Numerical solution of the system \f$ Ax=b \f$ by using Jacobi iteration.

The iteration is given by \f$ x^{(k+1)}=(I-Q^{-1}A)x^{(k)}+Q^{-1}b \f$ where
Q is the diagonals of A.*/
long int    CHI_MATH::SolveUsingJacobiIteration(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b, double tolerance, bool initialZero)
{
	long int n = A->rowCount;
	//===================================================== Defining Iteration vectors
	CHI_MATH_MATRIX* x_kp1 = new CHI_MATH_MATRIX(A->rowCount, 1);
	CHI_MATH_MATRIX* x_k = new CHI_MATH_MATRIX(A->rowCount, 1);
	if (!initialZero)
	{
		for (long int i = 1; i <= A->rowCount; i++)
		{
			x_k->sij(i,1,b->ij(i));
		}
	}

	//===================================================== Precalculating some values
	//CHI_MATH_MATRIX* mult1 = new CHI_MATH_MATRIX;
	//CHI_MATH_MATRIX* subt1 = new CHI_MATH_MATRIX;

	//===================================================== Calculating initial residual
	CHI_MATH_MATRIX* mult2 = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* subt2 = new CHI_MATH_MATRIX;
	mult2->OpMult(A, x_k);
	subt2->OpSubtract(b, mult2);
	double residual0 = subt2->Norm2();
	double residualk = 1.0;


	//===================================================== Iterating
	double sum1;
	double sum2;
	long int m = 1;
	bool firstIteration = true;
	while (((residualk / residual0)>tolerance) || (firstIteration))
	{
		if (firstIteration){ firstIteration = false; }
		else { x_k->Equals(x_kp1); }

		for (long int i = 1; i <= n; i++)
		{
			sum1 = 0.0;
			for (long int j = 1; j <= (i-1); j++)
			{
				sum1 += A->ij(i, j)*x_k->ij(j);
			}
			sum2 = 0.0;
			for (long int j = (i+1); j <= n; j++)
			{
				sum2 += A->ij(i, j)*x_k->ij(j);
			}
			x_kp1->sij(i,1, (b->ij(i) - sum1 - sum2) / A->ij(i, i));
		}

		//======================================= Calculating Residual k
		mult2->OpMult(A, x_kp1);
		subt2->OpSubtract(b, mult2);
		residualk = subt2->Norm2();

		m++;
	}

	x = x_kp1;

	return m;
}



//############################################################################# SOR iteration solver
/** This routine solves the system of equations \f$ Ax=b \f$ using the 
Successive Over Relaxation (SOR) iteration method.*/
long int CHI_MATH::SolveUsingSORIteration(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b, double omega, double tolerance, bool initialZero)
{
	long int n = A->rowCount;

	//===================================================== Defining Iteration vectors
	CHI_MATH_MATRIX* x_kp1 = new CHI_MATH_MATRIX(A->rowCount, 1);
	CHI_MATH_MATRIX* x_k = new CHI_MATH_MATRIX(A->rowCount, 1);
	if (!initialZero)
	{
		for (long int i = 1; i <= A->rowCount; i++)
		{
			x_k->sij(i,1, b->ij(i));
		}
	}

	//===================================================== Performing iteration
	//CHI_MATH_MATRIX* subt1 = new CHI_MATH_MATRIX;

	//===================================================== Calculating initial residual
	CHI_MATH_MATRIX* mult2 = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* subt2 = new CHI_MATH_MATRIX;
	mult2->OpMult(A, x_k);
	subt2->OpSubtract(b, mult2);
	double residual0 = subt2->Norm2();
	double residualk = 1.0;

	//===================================================== Iterating
	double sum1;
	double sum2;
	long int m = 1;
	bool firstIteration = true;
	while (((residualk / residual0) > tolerance) || (firstIteration))
	{
		if (firstIteration){ firstIteration = false; }
		else { x_k->Equals(x_kp1); }

		for (long int i = 1; i <= n; i++)
		{
			sum1 = 0.0;
			sum2 = 0.0;
			for (long int j = 1; j <= (i - 1); j++)	
			{ 
				sum1 += A->ij(i, j)*x_kp1->ij(j); 
			}
			for (long int j = (i+1); j <= n; j++)		
			{ 
				sum2 += A->ij(i, j)*x_k->ij(j);	
			}

			x_kp1->sij(i,1,(1.0-omega)*x_k->ij(i) + (omega / A->ij(i, i))*(b->ij(i) - sum1 - sum2));
		}

		//======================================= Calculating Residual k
		mult2->OpMult(A, x_kp1);
		subt2->OpSubtract(b, mult2);
		residualk = subt2->Norm2();

		m++;
	}

	x = x_kp1;

	return m;
}




//############################################################################# Symmetric SOR Solver
/** Implementation of the Symmetric SOR solver.*/
long int CHI_MATH::SolveUsingSymSORIteration(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b, double omega, double tolerance, bool initialZero)
{
	long int n = A->rowCount;

	//===================================================== Defining Iteration vectors
	CHI_MATH_MATRIX* x_kp5 = new CHI_MATH_MATRIX(A->rowCount, 1);
	CHI_MATH_MATRIX* x_kp1 = new CHI_MATH_MATRIX(A->rowCount, 1);
	CHI_MATH_MATRIX* x_k = new CHI_MATH_MATRIX(A->rowCount, 1);
	if (!initialZero)
	{
		for (long int i = 1; i <= A->rowCount; i++)
		{
			x_k->sij(i,1, b->ij(i));
		}
	}

	//===================================================== Calculating initial residual
	CHI_MATH_MATRIX* mult2 = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* subt2 = new CHI_MATH_MATRIX;
	mult2->OpMult(A, x_k);
	subt2->OpSubtract(b, mult2);
	double residual0 = subt2->Norm2();
	double residualk = 1.0;

	//===================================================== Iterating
	double sum1;
	double sum2;
	long int m = 1;
	bool firstIteration = true;
	while (((residualk / residual0) > tolerance) || (firstIteration))
	{
		if (firstIteration){ firstIteration = false; }
		else { x_k->Equals(x_kp1); }

		//============================= Forward
		for (long int i = 1; i <= n; i++)
		{
			sum1 = 0.0;
			for (long int j = 1; j <= (i - 1); j++)
			{
				sum1 += A->ij(i, j)*x_kp5->ij(j);
			}

			sum2 = 0.0;
			for (long int j = (i + 1); j <= n; j++)
			{
				sum2 += A->ij(i, j)*x_k->ij(j);
			}
			x_kp5->sij(i,1,(1.0 - omega)*x_k->ij(i) + (omega / A->ij(i, i))*(b->ij(i) - sum1 - sum2));
		}
		//============================= Backward
		for (long int i = n; i >= 1; i--)
		{
			sum1 = 0.0;
			for (long int j = 1; j <= (i - 1); j++)
			{
				sum1 += A->ij(i, j)*x_kp5->ij(j);
			}
			sum2 = 0.0;
			for (long int j = (i + 1); j <= n; j++)
			{
				sum2 += A->ij(i, j)*x_kp1->ij(j);
			}
			x_kp1->sij(i,1, (1.0 - omega)*x_kp5->ij(i) + (omega / A->ij(i, i))*(b->ij(i) - sum1 - sum2));
		}

		//======================================= Calculating Residual k
		mult2->OpMult(A, x_kp1);
		subt2->OpSubtract(b, mult2);
		residualk = subt2->Norm2();

		m++;
	}

	x = x_kp1;

	return m;
}

//############################################################################# Solve using steepest gradient
/** Solve the system using steepest gradient method.*/
long int CHI_MATH::SolveUsingSteepestDescent(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b, double tolerance, bool initialZero)
{
	long int n = A->rowCount;

	//===================================================== Defining Iteration vectors
	CHI_MATH_MATRIX* x_kp1 = new CHI_MATH_MATRIX(A->rowCount, 1);
	CHI_MATH_MATRIX* x_k = new CHI_MATH_MATRIX(A->rowCount, 1);
	if (!initialZero)
	{
		for (long int i = 1; i <= A->rowCount; i++)
		{
			x_k->sij(i,1, b->ij(i));
		}
	}

	//===================================================== Calculating initial residual
	CHI_MATH_MATRIX* Axk = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* v = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* Av = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* tv = new CHI_MATH_MATRIX;

	CHI_MATH_MATRIX* mult2 = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* subt2 = new CHI_MATH_MATRIX;
	mult2->OpMult(A, x_k);
	subt2->OpSubtract(b, mult2);
	double residual0 = subt2->Norm2();
	double residualk = 1.0;

	long int m = 1;
	bool firstIteration = true;
	while (((residualk / residual0) > tolerance) || (firstIteration))
	{
		if (firstIteration){ firstIteration = false; }
		else { x_k->Equals(x_kp1); }

		//======================================= Calculate direction vector
		Axk->OpMult(A, x_k);
		v->OpSubtract(b, Axk);
		Av->OpMult(A, v);

		//======================================= Calculating inner product
		double sum1 = 0.0;
		double sum2 = 0.0;
		for (long int k = 1; k <= n; k++)
		{
			sum1 += v->ij(k)*v->ij(k);
			sum2 += v->ij(k)*Av->ij(k);
		}
		double t = sum1 / sum2;

		//======================================= Calculating xkp1
		tv->OpMult(v, t);
		x_kp1->OpAdd(x_k, tv);

		//======================================= Calculating Residual k
		mult2->OpMult(A, x_kp1);
		subt2->OpSubtract(b, mult2);
		residualk = subt2->Norm2();

		m++;
	}

	x = x_kp1;

	return m;
}






//############################################################################# Solve using steepest gradient
/** Solve the system using steepest gradient method.*/
long int CHI_MATH::SolveUsingConjugantGrad(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b, double tolerance, bool initialZero)
{
	//long int n = A->rowCount;

	//===================================================== Defining Iteration vectors
	CHI_MATH_MATRIX* x_kp1 = new CHI_MATH_MATRIX(A->rowCount, 1);
	CHI_MATH_MATRIX* x_k = new CHI_MATH_MATRIX(A->rowCount, 1);
	if (!initialZero)
	{
		for (long int i = 1; i <= A->rowCount; i++)
		{
			x_k->sij(i,1,b->ij(i));
		}
	}

	//===================================================== Initializing storage
	CHI_MATH_MATRIX* v = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* z = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* Av = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* tv = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* tz = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* rtemp = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* dcv = new CHI_MATH_MATRIX;

	//===================================================== Calculating initial residual
	CHI_MATH_MATRIX* rmult = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* r = new CHI_MATH_MATRIX;

	CHI_MATH_MATRIX* rmult2 = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* r2 = new CHI_MATH_MATRIX;

	rmult->OpMult(A, x_k);
	r->OpSubtract(b, rmult);

	double residual0 = r->Norm2();
	double residualk = 1.0;

	v = r;

	double c = 0.0;	for (long int k = 1; k <= v->rowCount; k++){ c += r->ij(k)*r->ij(k); }
	double d = 0.0;
	double t = 0.0;

	long int m = 1;
	bool firstIteration = true;
	while (((residualk / residual0) > tolerance) || (firstIteration))
	{
		if (firstIteration){ firstIteration = false; }
		else { x_k->Equals(x_kp1); }

		Av->OpMult(A, v);
		z = Av;

		double sum_cz = 0.0; for (long int k = 1; k <= z->rowCount; k++){ sum_cz += v->ij(k)*z->ij(k); }
		t = c / sum_cz;

		tv->OpMult(v, t);

		x_kp1->OpAdd(x_k, tv);

		tz->OpMult(z, t);

		rtemp->OpSubtract(r, tz);
		r = rtemp;

		d = 0.0; for (long int k = 1; k <= v->rowCount; k++){ d += r->ij(k)*r->ij(k); }

		dcv->OpMult(v, d / c);

		v->OpAdd(r, dcv);

		c = d;

		//============================ Calculating residual
		rmult2->OpMult(A, x_kp1);
		r2->OpSubtract(b, rmult2);

		residualk = r2->Norm2();

		m++;
	}

	//x = x_kp1;
    for (long int k = 1; k <= x_kp1->rowCount; k++)
    {
        x->sij(k,1,x_kp1->ij(k,1));
    }
    
    //================================= Clearing Vectors        
    
    x_kp1  ->DestroyEntries();
    x_k    ->DestroyEntries();
    v      ->DestroyEntries();
    z      ->DestroyEntries();
    Av     ->DestroyEntries();
    tv     ->DestroyEntries();
    tz     ->DestroyEntries();
    rtemp  ->DestroyEntries();
    dcv    ->DestroyEntries();
    rmult  ->DestroyEntries();
    r      ->DestroyEntries();
    rmult2 ->DestroyEntries();
    r2     ->DestroyEntries();
    
    
    
    
    
    
    
    

	return m;
}








//############################################################################# Solve using conjugate gradient
/** Solve the system using steepest gradient method.*/
long int CHI_MATH::SolveSparseUsingConjugantGrad(CHI_MATH_SPARSE_MATRIX* A, CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b, double tolerance, bool initialZero)
{
	//long int n = A->rowCount;

	//===================================================== Defining Iteration vectors
	CHI_MATH_MATRIX* x_kp1 = new CHI_MATH_MATRIX(A->rowCount, 1);
	CHI_MATH_MATRIX* x_k = new CHI_MATH_MATRIX(A->rowCount, 1);
	if (!initialZero)
	{
		for (long int i = 1; i <= A->rowCount; i++)
		{
			x_k->sij(i,1, b->ij(i));
		}
	}

	//===================================================== Initializing storage
	CHI_MATH_MATRIX* v = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* z = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* Av = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* tv = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* tz = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* rtemp = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* dcv = new CHI_MATH_MATRIX;

	//===================================================== Calculating initial residual
	CHI_MATH_MATRIX* rmult = new CHI_MATH_MATRIX;
	CHI_MATH_MATRIX* r = new CHI_MATH_MATRIX;

	//CHI_MATH_MATRIX* rmult2 = new CHI_MATH_MATRIX;
	//CHI_MATH_MATRIX* r2 = new CHI_MATH_MATRIX;

	rmult->OpMult(A, x_k);
	r->OpSubtract(b, rmult);

	double residual0 = r->Norm2();
	double residualk = 1.0;

	v = r;

	double c = 0.0;	for (long int k = 1; k <= v->rowCount; k++){ c += r->ij(k)*r->ij(k); }
	double d = 0.0;
	double t = 0.0;

	long int m = 1;
	bool firstIteration = true;
	while (((residualk / residual0) > tolerance) || (firstIteration))
	{
		if (firstIteration){ firstIteration = false; }
		else { x_k->Equals(x_kp1); }

		Av->OpMult(A, v);
		z = Av;

		double sum_cz = 0.0; for (long int k = 1; k <= z->rowCount; k++){ sum_cz += v->ij(k)*z->ij(k); }
		t = c / sum_cz;

		tv->OpMult(v, t);

		x_kp1->OpAdd(x_k, tv);

		tz->OpMult(z, t);

		rtemp->OpSubtract(r, tz);
		r = rtemp;

		d = 0.0; for (long int k = 1; k <= v->rowCount; k++){ d += pow(r->ij(k),2); }

		dcv->OpMult(v, d / c);

		v->OpAdd(r, dcv);

		c = d;

		//============================ Calculating residual
		//rmult2->OpMult(A, x_kp1);
		//r2->OpSubtract(b, rmult2);

		//residualk = r2->Norm2();
		residualk = r->Norm2();

		m++;
	}
	
	x = x_kp1;

	return m;
}













