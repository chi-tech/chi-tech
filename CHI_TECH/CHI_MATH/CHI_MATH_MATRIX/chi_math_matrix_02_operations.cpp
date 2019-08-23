#include"chi_math_matrix.h"



//############################################################################# Matrix multiplication
/** This operator function handles matrix multiplications.

The basic formula for matrix multiplication, for \f$ A \in \mathbf{R}^{m\times p} \f$
and \f$ B \in \mathbf{R}^{p\times n} \f$ the multiplication is given by:
\f[
AB \in \mathbf{R}^{m\times n} \to
(AB)_{ij} = \sum_{k=1}^p a_{ik} b_{kj} \ \ (1\le i \le m,1\le j \le n)
\f]
*/
CHI_MATH_MATRIX* CHI_MATH_MATRIX::Mult(CHI_MATH_MATRIX* B)
{
	long int m;
	long int n;
	long int p;
	CHI_MATH_MATRIX* A = this;
	CHI_MATH_MATRIX* AB;

	//===================================================== Checking size compatibility
	if (A->colCount != B->rowCount)
	{
		//MessageBoxA(NULL, "Invalid matrix indexes for mult", "error", MB_OK);
		return 0;
	}
	else
	{
		m = A->rowCount;
		p = A->colCount;
		n = B->colCount;
		AB = new CHI_MATH_MATRIX(m, n);
		this->garbage.PushItem(AB);
	}

	//===================================================== Multiplying
	for (long int i = 1; i <= m; i++)
	{
		for (long int j = 1; j <= n; j++)
		{
			
			for (long int k = 1; k <= p; k++)
			{
				AB->sij(i, j, AB->ij(i, j)+ A->ij(i, k)*B->ij(k, j));
			}
		}
	}

	return AB;
}



//############################################################################# Constant multiplication
/** This function handles a matrix multiplied with a constant.*/
CHI_MATH_MATRIX* CHI_MATH_MATRIX::Mult(double alpha)
{
	CHI_MATH_MATRIX* alphaA = new CHI_MATH_MATRIX(this->rowCount, this->colCount);
	this->garbage.PushItem(alphaA);
	//===================================================== Multiplying
	for (long int i = 1; i <= this->rowCount; i++)
	{
		for (long int j = 1; j <= this->colCount; j++)
		{
			alphaA->sij(i, j,this->ij(i, j)*alpha);
		}
	}

	return alphaA;
}


//############################################################################# Matrix addition
/** This function handles a matrix added to another.*/
CHI_MATH_MATRIX* CHI_MATH_MATRIX::Add(CHI_MATH_MATRIX* B)
{
	//===================================================== Checking size compatibility
	if ((this->colCount != B->colCount) || (this->rowCount != B->rowCount))
	{
		//MessageBoxA(NULL, "Incompatible matrix sizes for A and B", "error", MB_OK);
		return 0;
	}

	CHI_MATH_MATRIX* ApB = new CHI_MATH_MATRIX(this->rowCount, this->colCount);
	this->garbage.PushItem(ApB);
	//===================================================== Adding
	for (long int i = 1; i <= this->rowCount; i++)
	{
		for (long int j = 1; j <= this->colCount; j++)
		{
			ApB->sij(i, j,this->ij(i, j) + B->ij(i, j));
		}
	}

	return ApB;
}

//############################################################################# Matrix subtraction
/** This function handles a matrix subtracted from another.*/
CHI_MATH_MATRIX* CHI_MATH_MATRIX::Sub(CHI_MATH_MATRIX* B)
{
	//===================================================== Checking size compatibility
	if ((this->colCount != B->colCount) || (this->rowCount != B->rowCount))
	{
		//MessageBoxA(NULL, "Incompatible matrix sizes for A and B", "error", MB_OK);
		return 0;
	}

	CHI_MATH_MATRIX* AmB = new CHI_MATH_MATRIX(this->rowCount, this->colCount);
	this->garbage.PushItem(AmB);
	//===================================================== Subtracting
	for (long int i = 1; i <= this->rowCount; i++)
	{
		for (long int j = 1; j <= this->colCount; j++)
		{
			AmB->sij(i, j,this->ij(i, j) - B->ij(i, j));
		}
	}

	return AmB;
}



//############################################################################# Transpose
/** Returns the transpose of the matrix.*/
CHI_MATH_MATRIX* CHI_MATH_MATRIX::Transpose()
{
	if (!this->transposeCalculated)
	{
		this->transpose = new CHI_MATH_MATRIX(this->colCount, this->rowCount);
		this->transposeCalculated = true;
		for (long int i = 1; i <= this->rowCount; i++)
		{
			for (long int j = 1; j <= this->colCount; j++)
			{
				this->transpose->sij(j, i,this->ij(i, j));
			}
		}
		
		return this->transpose;
	}
	else
	{
		for (long int i = 1; i <= this->rowCount; i++)
		{
			for (long int j = 1; j <= this->colCount; j++)
			{
				this->transpose->sij(j, i,this->ij(i, j));
			}
		}
		return this->transpose;
	}
}



//############################################################################# Operator =
void CHI_MATH_MATRIX::Equals(CHI_MATH_MATRIX* B)
{
	if (!this->initialized)
	{
		this->SetSizeZero(B->rowCount, B->colCount);
	}
	else
	{
		if ((this->rowCount != B->rowCount) || (this->colCount != B->colCount))
		{
			this->Resize(B->rowCount, B->colCount);
		}
	}

	for (long int i = 1; i <= this->rowCount; i++)
	{
		for (long int j = 1; j <= this->colCount; j++)
		{
			this->sij(i, j, B->ij(i, j));
		}
	}
}



//############################################################################# Norm_2 calculation
/** This operation calculates the euclidian norm of a matrix. In the special case
where the matrix is a single column (i.e. a vector), the euclidian norm is defined
as:
\f[
	||x||=\sum_{i=1}^{n} {|{x_i}^2|}^{1/2}
\f] */
double CHI_MATH_MATRIX::Norm2(long int column)
{
	double sum = 0.0;
	for (long int k = 1; k <= this->rowCount; k++)
		{
		sum += this->ij(k, column)*this->ij(k, column);
		}
		sum = sqrt(sum);
	
	return sum;
}





//############################################################################# Determine Maximum
/** Determines the maximum entry in the matrix.*/
double CHI_MATH_MATRIX::Max()
{
	double maximum = this->ij(1, 1);
	for (long int i = 1; i <= this->rowCount; i++)
	{
		for (long int j = 1; j <= this->colCount; j++)
		{
			if (this->ij(i, j) > maximum)
			{
				maximum = this->ij(i, j);
			}
		}
	}
	return maximum;
}




//############################################################################# Determine Minimum
/** Determines the minimum entry in the matrix.*/
double CHI_MATH_MATRIX::Min()
{
	double minimum = this->ij(1, 1);
	for (long int i = 1; i <= this->rowCount; i++)
	{
		for (long int j = 1; j <= this->colCount; j++)
		{
			if (this->ij(i, j) < minimum)
			{
				minimum = this->ij(i, j);
			}
		}
	}
	return minimum;
}
