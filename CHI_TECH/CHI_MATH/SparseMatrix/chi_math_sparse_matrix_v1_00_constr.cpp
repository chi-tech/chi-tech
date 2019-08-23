#include"chi_math_sparse_matrix_v1.h"

//######################################################### Default Constructor
/** Defines a zero matrix of size rows by columns.*/
CHI_MATH_SPARSE_MATRIX::CHI_MATH_SPARSE_MATRIX()
{
	this->initialized = false;
}

//######################################################### Constructor
/** Creates an empty matrix.*/
CHI_MATH_SPARSE_MATRIX::CHI_MATH_SPARSE_MATRIX(long int rows)
{

	this->compressedRows.threadProtectionEnabled = false;
	this->rowValues.threadProtectionEnabled = false;
	for (long int k = 0; k <= rows; k++)
	{
		CHI_VECTOR<long int>* newCompressedRow = new CHI_VECTOR < long int > ;
		CHI_VECTOR<double>* newCompressedRowVal = new CHI_VECTOR < double>;
		newCompressedRow->threadProtectionEnabled = false;
		newCompressedRowVal->threadProtectionEnabled = false;
		this->compressedRows.PushItem(newCompressedRow);
		this->rowValues.PushItem(newCompressedRowVal);
	}
	
	this->rowCount = rows;
	this->colCount = rows;
	this->initialized = true;
}