#include"chi_math_sparse_solver_cg.h"


//############################################################################# Constructor
CHI_MATH_SPARSE_SOLVER_CG::CHI_MATH_SPARSE_SOLVER_CG()
{
	this->tolerance = 1.0e-12;
	this->initialZero = true;
	this->initialValue = 0.0;
}

//############################################################################# Set initial guess
void CHI_MATH_SPARSE_SOLVER_CG::InitialGuess(double value)
{
	if (value != this->initialValue)
	{
		this->initialZero = false;
		this->initialValue = value;
	}
}

//############################################################################# Set initial guess
void CHI_MATH_SPARSE_SOLVER_CG::SetTolerance(double value)
{
	if ((value < 0.1) && (value>1.0e-24))
	{
		this->tolerance = value;
	}	
}

//############################################################################# Compute initial stuff
/** This routine allocates all the resources needed for the solve.*/
void CHI_MATH_SPARSE_SOLVER_CG::Compute(CHI_MATH_SPARSE_MATRIX* Ain)
{
	//===================================================== Defining Iteration vectors
	this->x_kp1 = new CHI_MATH_MATRIX(Ain->rowCount, 1);
	this->x_k = new CHI_MATH_MATRIX(Ain->rowCount, 1);

	//===================================================== Initializing storage
	this->v     = new CHI_MATH_MATRIX;
	this->z     = new CHI_MATH_MATRIX;
	this->Av    = new CHI_MATH_MATRIX;
	this->tv    = new CHI_MATH_MATRIX;
	this->tz    = new CHI_MATH_MATRIX;
	this->rtemp = new CHI_MATH_MATRIX;
	this->dcv   = new CHI_MATH_MATRIX;

	//===================================================== Calculating initial residual
	this->rmult = new CHI_MATH_MATRIX;
	this->r     = new CHI_MATH_MATRIX;

	this->rmult2 = new CHI_MATH_MATRIX;
	this->r2 = new CHI_MATH_MATRIX;

	this->A = Ain;
}


//############################################################################# Solve
/** Solve the sparse system using the conjugate gradient solver.*/
int CHI_MATH_SPARSE_SOLVER_CG::Solve(CHI_MATH_MATRIX*& x, CHI_MATH_MATRIX* b)
{
	//long int n = A->rowCount;

	//===================================================== Defining Iteration vectors
	if (!initialZero)
	{
		for (long int i = 1; i <= A->rowCount; i++)
		{
			x_k->sij(i,1,b->ij(i));
		}
	}
	else
	{
		for (long int i = 1; i <= A->rowCount; i++)
		{
			x_k->sij(i,1, 0.0);
		}
	}
	
	//===================================================== Calculating initial residual
	rmult->OpMult(A, x_k);
	r->OpSubtract(b, rmult);

	double residual0 = r->Norm2();
	double residualk = 1.0;

	v->Equals(r);

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
		z->Equals(Av);

		double sum_cz = 0.0; for (long int k = 1; k <= z->rowCount; k++){ sum_cz += v->ij(k)*z->ij(k); }
		t = c / sum_cz;

		tv->OpMult(v, t);

		x_kp1->OpAdd(x_k, tv);

		tz->OpMult(z, t);

		rtemp->OpSubtract(r, tz);
		r->Equals(rtemp);

		d = 0.0; for (long int k = 1; k <= v->rowCount; k++){ d += pow(r->ij(k), 2); }

		dcv->OpMult(v, d / c);

		v->OpAdd(r, dcv);

		c = d;

		//============================ Calculating residual
		residualk = r->Norm2();

		m++;
	}

	

	x->Equals(x_kp1);
	/*for (long int i = 1; i <= x_kp1->rowCount; i++)
	{		
		x->sij(i, 1, x_kp1->ij(i, 1));
	}*/
	//x_kp1->WriteToFile("../ZZZ.txt");
	//x->WriteToFile("../ZZZ.txt");

	/*for (long int i = 1; i <= rtemp->rowCount; i++)
	{
		rtemp->sij(i, 1, 0.0);
		x_k->sij(i, 1, 0.0);
		x_kp1->sij(i, 1, 0.0);
		v		->sij(i, 1, 0.0);
		z		->sij(i, 1, 0.0);
		Av		->sij(i, 1, 0.0);
		tv		->sij(i, 1, 0.0);
		tz		->sij(i, 1, 0.0);
		rtemp	->sij(i, 1, 0.0);
		dcv		->sij(i, 1, 0.0);
		rmult	->sij(i, 1, 0.0);
		r		->sij(i, 1, 0.0);
	}*/
	//x->WriteToFile("../ZZZ.txt");
	return m;
}
