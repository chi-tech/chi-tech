#include "chi_math_solver_cg.h"



//############################################################################# Initialize
/** Initializes all the utility vectors.
\author Jan*/
void CHI_MATH_SOLVER_CG::Initialize(long int size)
{
  x_kp1 =new CHI_MATH_MATRIX(size,1); 
  x_k   =new CHI_MATH_MATRIX(size,1);   
  v     =new CHI_MATH_MATRIX(size,1);    
  z     =new CHI_MATH_MATRIX(size,1);     
  Av    =new CHI_MATH_MATRIX(size,1);    
  tv    =new CHI_MATH_MATRIX(size,1);    
  tz    =new CHI_MATH_MATRIX(size,1);    
  rtemp =new CHI_MATH_MATRIX(size,1); 
  dcv   =new CHI_MATH_MATRIX(size,1);   
  rmult =new CHI_MATH_MATRIX(size,1); 
  r     =new CHI_MATH_MATRIX(size,1);     
  rmult2=new CHI_MATH_MATRIX(size,1);
  r2    =new CHI_MATH_MATRIX(size,1);    
}




//############################################################################# Solve
/** Uses the conjugate gradient method to solve the system.
\author Jan*/
long int  CHI_MATH_SOLVER_CG::Solve(CHI_MATH_MATRIX* A, 
                                    CHI_MATH_MATRIX*& x, 
                                    CHI_MATH_MATRIX* b, 
                                    double tolerance, 
                                    bool initialZero)
{

	//===================================================== Defining Iteration vectors
	if (!initialZero)
	{
		for (long int i = 1; i <= A->rowCount; i++)
		{
			x_k->sij(i,1,b->ij(i));
		}
	}

	//===================================================== Calculating initial residual
	rmult->OpMult(A, x_k);
	r->OpSubtract(b, rmult);

	double residual0 = r->Norm2();
	double residualk = 1.0;

    for (long int k = 1; k <= v->rowCount; k++)
    {
        v->sij(k,1,r->ij(k,1));
    }
	//v = r;

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
		//z = Av;
        for (long int k = 1; k <= Av->rowCount; k++)
        {
            z->sij(k,1,Av->ij(k,1));
        }

		double sum_cz = 0.0; for (long int k = 1; k <= z->rowCount; k++){ sum_cz += v->ij(k)*z->ij(k); }
		t = c / sum_cz;

		tv->OpMult(v, t);

		x_kp1->OpAdd(x_k, tv);

		tz->OpMult(z, t);

		rtemp->OpSubtract(r, tz);
		//r = rtemp;
        for (long int k = 1; k <= rtemp->rowCount; k++)
        {
            r->sij(k,1,rtemp->ij(k,1));
        }

		d = 0.0; for (long int k = 1; k <= v->rowCount; k++){ d += r->ij(k)*r->ij(k); }

		dcv->OpMult(v, d / c);

		v->OpAdd(r, dcv);

		c = d;

		//============================ Calculating residual
		rmult2->OpMult(A, x_kp1);
		r2->OpSubtract(b, rmult2);

		residualk = r2->Norm2();

		m++;

		if (m>10000)
		{
			std::cout<<"CG solver did not converge\n";
			break;
		}
	}

    //================================ Exporting the values
    for (long int k = 1; k <= x_kp1->rowCount; k++)
    {
        x->sij(k,1,x_kp1->ij(k,1));
    }


	return m;
}
