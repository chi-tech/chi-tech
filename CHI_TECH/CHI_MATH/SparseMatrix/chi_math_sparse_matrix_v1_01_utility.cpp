#include <cstdlib>
#include"chi_math_sparse_matrix_v1.h"

//#include"../../chi_tech.h"


//############################################################################# Set size with zeros
/** Defines a zero matrix of size rows by columns.*/
void CHI_MATH_SPARSE_MATRIX::SetSizeZero(long int rows)
{
	if (this->initialized)
	{
		//MessageBox(NULL, "Matrix already initialized", "Error", MB_ICONERROR);
		return;
	}



	this->compressedRows.threadProtectionEnabled = false;
	this->rowValues.threadProtectionEnabled = false;
	for (long int k = 0; k <= rows; k++)
	{
		CHI_VECTOR<long int>* newCompressedRow = new CHI_VECTOR < long int >;
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


//############################################################################# Create from text file (spaces)
/** Creates a multi column matrix from a space seperated file.

\note This algorithm is memory leak safe, except for the creation and filling of the matrix.*/
bool CHI_MATH_SPARSE_MATRIX::CreateFromText(char* fileName)
{
	char               fileLine[2000];
	char*              word;
	//int                numberOfValues;
	double             rawValue1 = 0.0;

	//===================================================== Opening the file
	FILE* file;
	file = fopen(fileName, "r");
	if (file == NULL)
	{ 
		//MessageBox(NULL, "File could not be read.", fileName, MB_ICONERROR); 
		return false; 
	}

	//===================================================== Reading every line and determining size
	int x = 0;
	int y = 0;
	int xmax = 0;
	int ymax = 0;
	while (feof(file) == 0)
	{
		fgets(fileLine, 1999, file);    //Reads a single line
		if (feof(file) == 0)
		{
			word = strtok(fileLine, " ");
			if (word != NULL)  {x = 1; y++;}
			while ((word != NULL) && (word[0] != char(10)))
			{
				sscanf(word, "%lf\n", &rawValue1);
				x++;
				word = strtok(NULL, " ");
			}
			x--;
			if (x > xmax) xmax = x;
			if (y > ymax) ymax = y;
		}
	}
	rewind(file);

	//===================================================== Creates the matrix from the file
	this->SetSizeZero(ymax);
	x = 0;
	y = 0;
	while (feof(file) == 0)
	{
		fgets(fileLine, 1999, file);    //Reads a single line
		if (feof(file) == 0)
		{
			word = strtok(fileLine, " ");
			if (word != NULL)  {x = 1; y++;}
			while ((word != NULL) && (word[0] != char(10)))
			{
				sscanf(word, "%lf\n", &rawValue1);
				if (fabs(rawValue1) > 1.0e-50)
				{
					this->s_ij(y, x, rawValue1);
				}
				x++;
				word = strtok(NULL, " ");
			}
			x--;
		}
	}
	fclose(file);

	return true;
}




//############################################################################# Setting a value
/** Creates and assigns a value to an entry. */
void CHI_MATH_SPARSE_MATRIX::s_ij(long int row, long int col,double value)
{
	if ((row>this->rowCount) || (col>this->colCount) || (row<1) || (col<1))
	{
		//MessageBox(NULL, "Invalid indexer", "Error", MB_ICONERROR);
		return;
	}
	long int* legalColumn;
	double*   legalValue;
    if (fabs(value) > 1.0e-50)
	{
		legalColumn = new long int; *legalColumn = col;
		legalValue = new double; *legalValue = value;
		this->compressedRows.GetItem(row - 1)->PushItem(legalColumn);
		this->rowValues.GetItem(row - 1)->PushItem(legalValue);
	}

}

//############################################################################# Accessing a value
/** The entries of a matrix are private and therefore,
this routine is used to gain access to the values. Will return
a reference to \f$ a_{1,1} \f$.*/
double CHI_MATH_SPARSE_MATRIX::ij(long int row, long int col)
{
	if ((row>this->rowCount) || (col>this->colCount) || (row<1) || (col<1))
	{
		//MessageBox(NULL, "Invalid indexer", "Error", MB_ICONERROR);
		return 0.0;
	}

	CHI_VECTOR<long int>* rowCol = this->compressedRows.GetItem(row - 1);
	CHI_VECTOR<double>* rowVal = this->rowValues.GetItem(row - 1);

	for (int k = 0; k < (rowCol->itemCount - 1); k++)
	{
		if (col == *rowCol->GetItem(k))
		{
			return *rowVal->GetItem(k);
		}
	}
	
	return 0.0;
}




//############################################################################# Compressing accessors
/** Finds all the legal columns for each row in order to optimize matrix multiplication.*/
void CHI_MATH_SPARSE_MATRIX::CompressRowWise()
{
	this->compressedRows.threadProtectionEnabled = false;
	for (long int i = 1; i <= this->rowCount;i++)
	{
		CHI_VECTOR<long int>* compressedRow = new CHI_VECTOR<long int>;
		compressedRow->threadProtectionEnabled = false;
		for (long int j = 1; j <= this->colCount; j++)
		{
			if (this->entries[i, j] != NULL)
			{
				long int* legalCol = new long int;
				*legalCol = j;
				compressedRow->PushItem(legalCol);
			}
		}
		this->compressedRows.PushItem(compressedRow);
	}
}
