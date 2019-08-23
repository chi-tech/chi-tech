#include"chi_math_matrix.h"



//############################################################################# Write matrix to file
/** This routine writes the matrix contents to a file.*/
void CHI_MATH_MATRIX::WriteToFile(char* fileName,int precision)
{
	std::ofstream file;
	file.open(fileName);
	//file.setf(std::ios::scientific);
	//file.precision(precision);

//	file << "[";
//	for (long int i = 1; i <= this->rowCount; i++)
//	{
//		for (long int j = 1; j <= this->colCount; j++)
//		{
//			if (this->ij(i, j) >= 0.0) { file << " "; }
//			file << this->ij(i, j) << " ";
//		}
//
//		file << " ;";
//		if (i!=this->rowCount) file << std::endl;
//	}
//	file << "]\n";
	file.close();
}



//############################################################################# Clean garbage vector
/** Clears the garbage vector.*/
void CHI_MATH_MATRIX::Clean()
{
	this->garbage.ClearVector();
}

//############################################################################# Destroys entries
/** Destroys internal entries data.*/
void CHI_MATH_MATRIX::DestroyEntries()
{
	this->garbage.ClearVector();
	for (long int i=0; i <= this->rowCount; i++)
	{
		delete[] this->entries[i];
	}
}


//############################################################################# Resizing
/**This routine resizes a matrix.
\warning Operation intensive.
*/
void CHI_MATH_MATRIX::Resize(long int rows, long int columns)
{
	this->DestroyEntries();
	this->entries = new double*[rows + 1];
	for (long int k = 0; k <= rows; k++)
	{
		this->entries[k] = new double[columns + 1];
		for (long int q = 0; q <= columns; q++)
		{
			this->entries[k][q] = 0.0;
		}
	}
	this->rowCount = rows;
	this->colCount = columns;
	this->initialized = true;
}
