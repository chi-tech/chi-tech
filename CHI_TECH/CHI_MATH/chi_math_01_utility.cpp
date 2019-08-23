#include "chi_math.h"



//############################################################################# Linear interpolation
/** Makes a linear interpolation given a lookup value and a lookup table.

This routine uses recursive subdivision to find the interpolation indices. The number 
of subdivisions that result in the most optimum interpolation behavior has been
determined using a progressive study of intervals.
- Subs = 100, 4286ms
- Subs =  10, 895ms
- Subs =   5, 618ms
- Subs =   4, 534ms
- Subs =   3, 506ms
- Subs =   2, 500ms

What is the math behind this? For n amount of values, with no recursive subdivision there would
be \f$2(n-1)\approx 2n\f$ comparison operations. Now, with each interval, we divide the interval into
\f$x\f$ amount of sub-intervals which results in \f$y\approx 2n/x\f$ sub-evaluations. Now suppose this 
happens recursively up to the point where there is less than or equal to \f$x\f$ lookup values left, we
will have the final linear lookup requiring \f$ x\f$ amount of lookups. Therefore, with level \f$l=1,\cdots,L\f$,
we get:
\f[
	y_{total}=\sum_{l=1}^L \biggr( \frac{2n}{x^l} \biggr) + x
\f]
So, from this formula, we can see that the minimum efficiency is linearly determined by \f$x\f$ and 
inversely by \f$n\f$. However, another competing effect is preset. This formulation does not account 
for the effect of randomized interval hits. In other words, the probality of hitting the appropriate 
interval in each substep is \f$\frac{1}{x}\f$ and in terms of this relation, the efficiency is 
inversely determined by \f$x\f$. From numerical analysis, for \f$n>1000\f$ it can be shown that the 
probability effect is a much more accurate reflection of the behavior, especially taking into 
consideration that successive multiplication of the probabilities results in a small number.

\warning Table indexColumn must be sorted in increasing order.*/
double CHI_MATH::LinearInterpolate(double lookUpValue, CHI_MATH_MATRIX* lookUpMatrix, int lookUpColumn, int indexColumn)
{
	//===================================================== If lookupValue <= smallest index
	if (lookUpValue <= lookUpMatrix->ij(1, indexColumn))
	{
		return lookUpMatrix->ij(1, lookUpColumn);
	}

	//===================================================== If lookupValue >= largest index
	if (lookUpValue >= lookUpMatrix->ij(lookUpMatrix->rowCount, indexColumn))
	{
		return lookUpMatrix->ij(lookUpMatrix->rowCount, lookUpColumn);
	}

	//===================================================== Determining optimization parameters
	int numberOfSubs = 2;
	//int numberOfLvls = 1;
	long int kfrom = 1;
	long int kto = lookUpMatrix->rowCount;

	//===================================================== Starting index scan
	int lvl = 1;
	bool minResolutionReached = false;
//#pragma warning(disable:4244)
	long int jumpAmount = floor((double)(kto - kfrom) / numberOfSubs);
	if (jumpAmount <= numberOfSubs) { minResolutionReached = true; }
	while (!minResolutionReached)
	{
		long int jumpAmount = floor((double)(kto - kfrom) / numberOfSubs);
//#pragma warning(default:4244)
		if (jumpAmount <= numberOfSubs) { minResolutionReached = true; break; }

		bool found = false;
		for (int k = 1; k <= numberOfSubs; k++)
		{
			if (((kfrom + jumpAmount) > kto) || (k == numberOfSubs))
			{
				if ((lookUpValue >= lookUpMatrix->ij(kfrom, indexColumn)) && (lookUpValue < lookUpMatrix->ij(kto, indexColumn)))
				{kto = kto;	}
			}
			else
			{
				if ((lookUpValue >= lookUpMatrix->ij(kfrom, indexColumn)) && (lookUpValue < lookUpMatrix->ij(kfrom + jumpAmount, indexColumn)))
				{
					kto = kfrom + jumpAmount;	found = true;
				}
			}
			if (!found) kfrom += jumpAmount;
			break;
		}

		lvl++;
	}

	//===================================================== Final interpolation
	//===================================================== If lookupValue <= smallest index
	if (lookUpValue <= lookUpMatrix->ij(kfrom, indexColumn))
	{
		return lookUpMatrix->ij(kfrom, lookUpColumn);
	}

	//===================================================== If lookupValue >= largest index
	if (lookUpValue >= lookUpMatrix->ij(kto, indexColumn))
	{
		return lookUpMatrix->ij(kto, lookUpColumn);
	}
	for (long int k = kfrom; k < kto; k++)
	{
		if ((lookUpValue >= lookUpMatrix->ij(k, indexColumn)) && (lookUpValue < lookUpMatrix->ij(k + 1, indexColumn)))
		{
	
			double returnValue;
			returnValue = (lookUpMatrix->ij(k + 1, lookUpColumn) - lookUpMatrix->ij(k, lookUpColumn));
			returnValue = returnValue / (lookUpMatrix->ij(k + 1, indexColumn) - lookUpMatrix->ij(k, indexColumn));
			returnValue = returnValue * (lookUpValue - lookUpMatrix->ij(k, indexColumn));
			returnValue = returnValue + lookUpMatrix->ij(k, lookUpColumn);
			return returnValue;
		}
	}

	

	//===================================================== Problem in lookup Value
	//MessageBox(NULL, "Failed linear interpolation", "Returning Zero", MB_ICONERROR);
	return lookUpMatrix->ij(lookUpMatrix->rowCount, lookUpColumn);
}








//############################################################################# Reindexes one table according to another
/** This routine index one matrix according to another. 

\warning Only to be used on nx2 tables.
*/
void CHI_MATH::ReIndexTable(CHI_MATH_MATRIX*& A, CHI_MATH_MATRIX* accordingToB)
{
	//===================================================== Create linInterp object for A
	/*CHI_MATH_LINEAR_INTERPOLATION linInterp;
	linInterp.Optimize(A);*/

	//===================================================== Create new matrix to replace A
	CHI_MATH_MATRIX* newM = new CHI_MATH_MATRIX(accordingToB->rowCount, accordingToB->colCount);
	
	//===================================================== Copy indexes from B and find corresponding values
	for (long int k = 1; k <= accordingToB->rowCount; k++)
	{
		newM->sij(k, 1,accordingToB->ij(k, 1));
		newM->sij(k, 2,this->LinearInterpolate(accordingToB->ij(k, 1), A, 2, 1));
	}
	A->DestroyEntries();
	A = newM;
}




//############################################################################# Create Test Matrix A
/** This routine creates a test matrix of a strictly diagonally dominant matrix. The matrix
has the following form:
\f[
A=
\begin{bmatrix}
    k_0+k_1  &    -k_1     & 0       & 0       & \cdots    & 0          & 0             \\
	   -k_1  &   k_1+k_2   & -k_2    & 0       & \cdots    & 0          & 0             \\
	 \vdots  &     \vdots  & \vdots  & \vdots  & \ddots    & \vdots     & \vdots         \\
	 0       &   0         & 0       & 0       & \cdots    & -k_{n-1}   & k_{n-1}+k_n
\end{bmatrix}
\ b=
\begin{bmatrix}
   0 \\
   0 \\
   \vdots \\
   k_n
\end{bmatrix}
\f]
In this form, the \f$k_n\f$ coefficients can either be formed from a linear representation (linear=true)
or a piecewise constant function. For the linear case: where \f$ h=\frac{1}{size+1}\f$ the values 
\f$ k_i = 1+(i+0.5)h \f$ are assigned. For the non-linear case: \f$ t=(i+0.5)h \f$ and 
\f$ k(t) = 1 \f$ when \f$ 0\le t < 0.5 \f$ and \f$ k(t)=K \f$ when \f$ 0.5\le t < 1 \f$. For 
\f$ K=2,100, 1000 \f$.

*/
void    CHI_MATH::CreateTestMatrixA(CHI_MATH_MATRIX*& A, CHI_MATH_MATRIX*& b, long int size, bool linear, double K)
{
	A = new CHI_MATH_MATRIX(size, size);
	b = new CHI_MATH_MATRIX(size, 1);

	//===================================================== Computing coeficient elements
	double* k = new double[size + 1];
	double  h = 1.0 / (size + 1.0);
	if (linear)
	{
		for (long int i = 0; i <= size; i++)
		{
			k[i] = (double)(1 + (i + 0.5)*h);
		}
	}
	else
	{
		for (long int i = 0; i <= size; i++)
		{
			double t = (double)((i + 0.5)*h);
			if ((t >= 0.0) && (t < 0.5))
			{
				k[i] = 1.0;
			}

			if ((t >= 0.5) && (t < 1.0))
			{
				k[i] = K;
			}
			
		}
	}

	//===================================================== Filling elements
	for (long int i = 1; i <= size; i++)
	{
		A->sij(i,i,k[i - 1] + k[i]);
		if (i > 1)    { A->sij(i - 1,i,  -1.0*k[i - 1]); }
		if (i < size) { A->sij(i + 1,i, -1.0*k[i]); }

		if (i < size) { b->sij(i,1, 0.0); }
		if (i == size) { b->sij(i,1,k[i]); }
	}
}



//############################################################################# Create Test Matrix B
/** This matrix uses lexicographical representation of the equation (\f$ h=\frac{1}{size+1}\f$):
\f[
	(4+h^2)x_{i,j} - x_{i-1,j} - x_{i+1,j} - x_{i,j-1} - x_{i,j+1} =h^2
\f]
\f[
	x_{0,j}=x_{n+1,j}=x_{i,0}=x_{i,n+1}=0
\f] */
void CHI_MATH::CreateTestMatrixB(CHI_MATH_MATRIX*& A, CHI_MATH_MATRIX*& b, long int size)
{
	long int n = size;
	long int n2 = size*size;
	A = new CHI_MATH_MATRIX(n2, n2);
	b = new CHI_MATH_MATRIX(n2, 1);

	//===================================================== Building A matrix
	double h = 1.0 / (n + 1.0);
	int l = 0;								 //Row index of a
	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			//=============================== Incrementing row index
			l += 1;

			//=============================== Determining lexicographic indices
			int m1 = (i - 1)*n + j;          //Lexicographic index of x i   j 
			int m2 = (i - 2)*n + j;                                 //x i-1 j
			int m3 = (i - 0)*n + j;                                 //x i+1 j
			int m4 = (i - 1)*n + j - 1;                             //x i   j-1
			int m5 = (i - 1)*n + j + 1;                             //x i   j+1

			//=============================== Setting result value
			b->sij(l,1, h*h);

			//=============================== Setting certain coefficients
			A->sij(l,m1, (4.0 + h*h));          //Setting coefficient of x i   j 

			if (i > 1) { A->sij(l,m2, -1.0); }  //Setting coefficient of x i-1 j 
			if (i < n) { A->sij(l,m3, -1.0); }  //Setting coefficient of x i+1 j 
			if (j > 1) { A->sij(l,m4, -1.0); }  //Setting coefficient of x i   j-1 
			if (j < n) { A->sij(l,m5, -1.0); }  //Setting coefficient of x i   j+1 
		}
	}
}



//############################################################################# Create test matrix C
/**Creates a random diagonally dominant matrix. Diagonals \f$a_{i,i}=i\f$ and diagonal elements are
\f$a_{i,j} = a_{i,i}/10^{|i-j|}\f$. \f$b_i=1\f$.*/
void  CHI_MATH::CreateTestMatrixC(CHI_MATH_MATRIX*& A, CHI_MATH_MATRIX*& b, long int size)
{
	long int n = size;
	A = new CHI_MATH_MATRIX(n, n);
	b = new CHI_MATH_MATRIX(n, 1);

	for (long int i = 1; i <= n; i++)
	{
		A->sij( i, i , (double)i);
		b->sij(i,1, 1.0);
		for (long int j = 1; j <= n; j++)
		{
			if (j != i)
			{
				A->sij(i, j, A->ij(i, i) / pow(10.0, (double)std::abs(i - j)));
			}
		}
	}
}
