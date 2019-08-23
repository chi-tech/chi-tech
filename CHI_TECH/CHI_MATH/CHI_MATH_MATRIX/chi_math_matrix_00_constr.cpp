#include"chi_math_matrix.h"

//######################################################### Default Constructor
/** Defines a zero matrix of size rows by columns.*/
CHI_MATH_MATRIX::CHI_MATH_MATRIX()
{
	this->initialized         = false;
	this->transposeCalculated = false;
}

//######################################################### Constructor
/** Defines a zero matrix of size rows by columns.*/
CHI_MATH_MATRIX::CHI_MATH_MATRIX(long int rows, long int columns)
{
	this->entries = new double*[rows + 1];
	for (long int k = 0; k <= rows; k++)
	{
		this->entries[k] = new double[columns + 1];
		for (long int n = 0; n <= columns; n++)
		{
			this->entries[k][n] = 0.0;
		}
	}
	this->rowCount = rows;
	this->colCount = columns;
	this->initialized = true;
}


