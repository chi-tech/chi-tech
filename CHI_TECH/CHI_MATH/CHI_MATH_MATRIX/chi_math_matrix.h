#ifndef CHI_MATH_MATRIX_H
#define CHI_MATH_MATRIX_H

#include "../chi_math_incdef.h"
#include "../../CHI_VECTOR/chi_vector.h"
#include "../SparseMatrix/chi_math_sparse_matrix_v1.h"

class CHI_MATH_MATRIX;
//template class CHI_VECTOR < CHI_MATH_MATRIX > ;
//######################################################### Class definition
/** Object for defining arbitrary sized matrices

A matrix can be created with a given size by using either a constructor,
or by setting its size later on.

Example:
\code
CHI_MATH_MATRIX Zt(2,3);   //Creates 2x3 matrix
CHI_MATH_MATRIX Zt;        //Creates empty matrix
Zt.SetSizeZero(2,3);       //Creates 2x3 matrix
\endcode

*/
class CHI_MATH_MATRIX
{
public:
	long int						rowCount;                  ///< Number of rows
	long int						colCount;                  ///< Number of columns
private:
	double**						entries;                   ///< Empty array of pointers
	bool							initialized;               ///< Flag indicating if it is initialized
	bool							transposeCalculated;       ///< Flag indicating if the transpose has been calculated
	CHI_MATH_MATRIX*				transpose;                 ///< Pre-calculated transpose vector
	CHI_VECTOR<CHI_MATH_MATRIX>     garbage;                   ///< Garbage vector
	CHI_VECTOR<CHI_MATH_MATRIX>     currentOperationsMatrix;   ///< A list of matrices that hold operations results
public:
	//00 Constr
						CHI_MATH_MATRIX();
						CHI_MATH_MATRIX(long int rows, long int columns);
	//01 Utility
	void				SetSizeZero(long int rows, long int columns);
	void				SetSizeValue(long int rows, long int columns, double value);
	bool				CreateFromText2Col(char* fileName, const char format[]="%lf,%lf");
	bool                CreateFromText(char* fileName);
	void                CreateFromLexicalVector(CHI_MATH_MATRIX* x, long int size);
	double				ij(long int row, long int col=1);
	void 				sij(long int row, long int col,const double value);

	//02 Operations
	CHI_MATH_MATRIX*    Mult(CHI_MATH_MATRIX* B);
	CHI_MATH_MATRIX*    Mult(double alpha);
	CHI_MATH_MATRIX*    Add(CHI_MATH_MATRIX* B);
	CHI_MATH_MATRIX*    Sub(CHI_MATH_MATRIX* B);
	CHI_MATH_MATRIX*    Transpose();
	CHI_MATH_MATRIX*    Inverse();
	void                Equals(CHI_MATH_MATRIX* B);
	double              Norm2(long int column=1);
	double              Max();
	double              Min();

	//03 Auxiliary
	void                WriteToFile(char* fileName,int precision=2);
	void                Clean();
	void                DestroyEntries();
	void                Resize(long int rows, long int columns);

	//04 OPoperations
	CHI_MATH_MATRIX*    OpMult(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* B);
	CHI_MATH_MATRIX*    OpMult(CHI_MATH_MATRIX* A, double alpha);
	CHI_MATH_MATRIX*    OpMult(CHI_MATH_SPARSE_MATRIX* A, CHI_MATH_MATRIX* B);
	CHI_MATH_MATRIX*    OpAdd(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* B);
	CHI_MATH_MATRIX*    OpSubtract(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* B, long int column);
	CHI_MATH_MATRIX*    OpSubtract(CHI_MATH_MATRIX* A, CHI_MATH_MATRIX* B);
};

#endif
