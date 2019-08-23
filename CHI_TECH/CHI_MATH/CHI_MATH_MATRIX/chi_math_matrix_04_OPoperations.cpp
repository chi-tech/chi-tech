#include"chi_math_matrix.h"	


//############################################################################# Matrix multiplication
/** Performs a self assigning multiplication. If the matrix dimension does not match
the required size, the matrix is reinitialized.
*/
CHI_MATH_MATRIX* CHI_MATH_MATRIX::OpMult(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* B)
{
	long int m;
	long int n;
	long int p;
	

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
		
		if (!this->initialized) 
		{
			this->SetSizeZero(m, n);
		}
		else
		{
			if ((this->rowCount!=m)||(this->colCount!=n))
			{
				this->Resize(m, n);
			}
		}
		
	}

	//===================================================== Multiplying
	for (long int i = 1; i <= m; i++)
	{
		for (long int j = 1; j <= n; j++)
		{
			this->sij(i, j, 0.0);
			for (long int k = 1; k <= p; k++)
			{
				this->sij(i, j, this->ij(i, j) + A->ij(i, k)*B->ij(k, j));
			}
		}
	}

	return this;
}




//############################################################################# Matrix multiplication
/** Performs a self assigning multiplication. If the matrix dimension does not match
the required size, the matrix is reinitialized.
*/
CHI_MATH_MATRIX* CHI_MATH_MATRIX::OpMult(CHI_MATH_MATRIX* A, double alpha)
{
	long int m;
	long int n;

	//===================================================== Checking size compatibility
	{
		m = A->rowCount;
		n = A->colCount;


		if (!this->initialized)
		{
			this->SetSizeZero(m, n);
		}
		else
		{
			if ((this->rowCount != m) || (this->colCount != n))
			{
				this->Resize(m, n);
			}
		}

	}

	//===================================================== Multiplying
	for (long int i = 1; i <= m; i++)
	{
		for (long int j = 1; j <= n; j++)
		{
			this->sij(i, j, A->ij(i, j)*alpha);
		}
	}

	return this;
}




//############################################################################# Matrix multiplication
/** Performs a self assigning multiplication. If the matrix dimension does not match
the required size, the matrix is reinitialized.
*/
CHI_MATH_MATRIX* CHI_MATH_MATRIX::OpMult(CHI_MATH_SPARSE_MATRIX* A, CHI_MATH_MATRIX* B)
{
	long int m;
	long int n;
	//long int p;


	//===================================================== Checking size compatibility
	if (A->colCount != B->rowCount)
	{
		//MessageBoxA(NULL, "Invalid matrix indexes for mult", "error", MB_OK);
		return 0;
	}
	else
	{
		m = A->rowCount;
		//p = A->colCount;
		n = B->colCount;

		if (!this->initialized)
		{
			this->SetSizeZero(m, n);
		}
		else
		{
			if ((this->rowCount != m) || (this->colCount != n))
			{
				this->Resize(m, n);
			}
		}

	}

	//===================================================== Multiplying
	long int legalCol = 0;
	for (long int i = 1; i <= m; i++)
	{
		for (long int j = 1; j <= n; j++)
		{
			this->sij(i, j,0.0);
			for (long int k = 0; k < A->compressedRows.GetItem(i-1)->itemCount; k++)
			{
				legalCol = *A->compressedRows.GetItem(i-1)->GetItem(k);
				this->sij(i, j, this->ij(i, j)+ *A->rowValues.GetItem(i - 1)->GetItem(k)*B->ij(legalCol, j));
			}
		}
	}

	return this;
}





//############################################################################# Matrix addition
/** This function handles a matrix added to another.*/
CHI_MATH_MATRIX* CHI_MATH_MATRIX::OpAdd(CHI_MATH_MATRIX* A,CHI_MATH_MATRIX* B)
{
	//===================================================== Checking size compatibility
	if ((A->colCount != B->colCount) || (A->rowCount != B->rowCount))
	{
		//MessageBoxA(NULL, "Incompatible matrix sizes for A and B", "error", MB_OK);
		return 0;
	}

	//===================================================== Checking Sizing
	long int m = A->rowCount;
	long int n = A->colCount;

	if (!this->initialized)
	{
		this->SetSizeZero(m, n);
	}
	else
	{
		if ((this->rowCount != m) || (this->colCount != n))
		{
			this->Resize(m, n);
		}
	}

	//===================================================== Adding
	for (long int i = 1; i <= this->rowCount; i++)
	{
		for (long int j = 1; j <= this->colCount; j++)
		{
			this->sij(i, j, A->ij(i, j) + B->ij(i, j));
		}
	}

	return this;
}


//############################################################################# Matrix sutraction
/** This function handles a matrix added to another.*/
CHI_MATH_MATRIX* CHI_MATH_MATRIX::OpSubtract(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* B, long int column)
{
	//===================================================== Checking size compatibility
	if ((A->colCount != B->colCount) || (A->rowCount != B->rowCount))
	{
		//MessageBoxA(NULL, "Incompatible matrix sizes for A and B", "error", MB_OK);
		return 0;
	}

	//===================================================== Checking Sizing
	long int m = A->rowCount;
	long int n = A->colCount;

	if (!this->initialized)
	{
		this->SetSizeZero(m, n);
	}
	else
	{
		if ((this->rowCount != m) || (this->colCount != n))
		{
			this->Resize(m, n);
		}
	}

	//===================================================== Subtracting
	for (long int i = 1; i <= this->rowCount; i++)
	{
		for (long int j = column; j <= column; j++)
		{
			this->sij(i, j, A->ij(i, j) - B->ij(i, j));
		}
	}

	return this;
}



//############################################################################# Matrix sutraction
/** This function handles a matrix added to another.*/
CHI_MATH_MATRIX* CHI_MATH_MATRIX::OpSubtract(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* B)
{
	//===================================================== Checking size compatibility
	if ((A->colCount != B->colCount) || (A->rowCount != B->rowCount))
	{
		//MessageBoxA(NULL, "Incompatible matrix sizes for A and B", "error", MB_OK);
		return 0;
	}

	//===================================================== Checking Sizing
	long int m = A->rowCount;
	long int n = A->colCount;

	if (!this->initialized)
	{
		this->SetSizeZero(m, n);
	}
	else
	{
		if ((this->rowCount != m) || (this->colCount != n))
		{
			this->Resize(m, n);
		}
	}

	//===================================================== Subtracting
	for (long int i = 1; i <= this->rowCount; i++)
	{
		for (long int j = 1; j <= this->colCount; j++)
		{
			this->sij(i, j,A->ij(i, j) - B->ij(i, j));
		}
	}

	return this;
}



