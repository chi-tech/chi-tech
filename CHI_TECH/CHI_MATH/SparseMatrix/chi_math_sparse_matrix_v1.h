#ifndef CHI_MATH_SPARSE_MATRIX_H
#define CHI_MATH_SPARSE_MATRIX_H

#include "../chi_math_incdef.h"
#include "../../CHI_VECTOR/chi_vector.h"

//template class CHI_VECTOR < long int > ;
//template class CHI_VECTOR < double >;
//template class CHI_VECTOR < CHI_VECTOR < long int > > ;
//template class CHI_VECTOR < CHI_VECTOR < double > >;

//######################################################### Class Definition
/** Sparse matrix object*/
class CHI_MATH_SPARSE_MATRIX
{
public:
	long int							rowCount;                  ///< Number of rows
	long int							colCount;                  ///< Number of columns
	CHI_VECTOR<CHI_VECTOR<long int> >    compressedRows;            ///< Storage of all legal columns
	CHI_VECTOR<CHI_VECTOR<double> >      rowValues;                 ///< Storage of all values

private:
	double***						entries;                   ///< Empty array of pointers
	bool							initialized;               ///< Flag indicating if it is initialized
	
public:
	//00 Constr
	CHI_MATH_SPARSE_MATRIX();
	CHI_MATH_SPARSE_MATRIX(long int rows);
	//01 Utility
	void				SetSizeZero(long int rows);
	bool                CreateFromText(char* fileName);
	void                s_ij(long int row, long int col, double value);
	double				ij(long int row, long int col);
	void                CompressRowWise();

	//02 Operations
	/*CHI_MATH_MATRIX*    Mult(CHI_MATH_MATRIX* B);
	CHI_MATH_MATRIX*    Mult(double alpha);
	CHI_MATH_MATRIX*    Add(CHI_MATH_MATRIX* B);
	CHI_MATH_MATRIX*    Sub(CHI_MATH_MATRIX* B);
	CHI_MATH_MATRIX*    Transpose();
	CHI_MATH_MATRIX*    Inverse();
	void                operator=(CHI_MATH_MATRIX* B);
	double              Norm2(long int column = 1);
	double              Max();
	double              Min();*/

	//03 Auxiliary
	/*void                WriteToFile(char* fileName, int precision = 2);
	void                Clean();
	void                DestroyEntries();
	void                Resize(long int rows, long int columns);*/

	//04 OPoperations
	/*CHI_MATH_MATRIX*    OpMult(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* B);
	CHI_MATH_MATRIX*    OpMult(CHI_MATH_MATRIX* A, double alpha);
	CHI_MATH_MATRIX*    OpAdd(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* B);
	CHI_MATH_MATRIX*    OpSubtract(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* B, long int column);
	CHI_MATH_MATRIX*    OpSubtract(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* B);*/
};


#endif
