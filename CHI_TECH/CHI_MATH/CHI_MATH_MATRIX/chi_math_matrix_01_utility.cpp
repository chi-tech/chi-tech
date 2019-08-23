#include"chi_math_matrix.h"


//############################################################################# Set size with zeros
/** Defines a zero matrix of size rows by columns.*/
void CHI_MATH_MATRIX::SetSizeZero(long int rows, long int columns)
{
	if (this->initialized)
	{
		//MessageBox(NULL, "Matrix already initialized", "Error", MB_ICONERROR);
		return;
	}

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
	this->initialized=true;
}

//############################################################################# Set size with values
/** Defines a matrix, of all given entries equal to value of size rows by columns.*/
void CHI_MATH_MATRIX::SetSizeValue(long int rows, long int columns, double value)
{
	if (this->initialized)
	{
		//MessageBox(NULL, "Matrix already initialized", "Error", MB_ICONERROR);
		return;
	}

	this->entries = new double*[rows + 1];
	for (long int k = 0; k <= rows; k++)
	{
		this->entries[k] = new double[columns + 1];
		for (long int n = 0; n <= columns; n++)
		{
			this->entries[k][n] = value;
		}
	}
	this->rowCount = rows;
	this->colCount = columns;
	this->initialized = true;
}




//############################################################################# Create from text file
/** Creates a 2 column matrix from a csv file.

\warning This algorithm creates a two vectors that are not cleared.

\todo Fix the uncleared vectors. They pop an error if you try.*/
bool CHI_MATH_MATRIX::CreateFromText2Col(char* fileName, const char* format)
{
	double             rawValue1=0.0;
	double             rawValue2=0.0;
	double*            pValue1;
	double*            pValue2;
	CHI_VECTOR<double> value1;
	CHI_VECTOR<double> value2;
	
	//===================================================== Opening the file
	FILE* file;
	file = fopen(fileName, "r");
	if (file == NULL){ //MessageBox(NULL, "File could not be read.", fileName, MB_ICONERROR); 
    return false; }
	
	//===================================================== Reading every line and filling arrays
	bool atEndofFile = false;
	int index = 0;
	while (!atEndofFile)
	{
		index = fscanf(file, format, &rawValue1, &rawValue2);
		if (index == EOF){ atEndofFile = true; }
		else
		{
			pValue1 = new double; *pValue1 = rawValue1;
			pValue2 = new double; *pValue2 = rawValue2;
			value1.PushItem(pValue1);
			value2.PushItem(pValue2);
		}
	}
	fclose(file);

	//===================================================== Creating the matrix
	this->SetSizeZero(value1.itemCount, 2);
	for (long int k = 0; k < value1.itemCount; k++)
	{
		this->entries[k + 1][1] = *value1.GetItem(k);
		this->entries[k + 1][2] = *value2.GetItem(k);
	}

	//===================================================== Cleaning the vectors
	//value1.ClearVector();
	//value2.ClearVector();

	return true;
}

//############################################################################# Create from text file (spaces)
/** Creates a multi column matrix from a space seperated file.

\note This algorithm is memory leak safe, except for the creation and filling of the matrix.*/
bool CHI_MATH_MATRIX::CreateFromText(char* fileName)
{
	char               fileLine[20000];
	char*              word;
	double             rawValue1 = 0.0;
	//double             rawValue2 = 0.0;
	CHI_VECTOR<double> value1;
	CHI_VECTOR<double> value2;

	//===================================================== Opening the file
	FILE* file;
	file = fopen(fileName, "r");
	if (file == NULL){ //MessageBox(NULL, "File could not be read.", fileName, MB_ICONERROR); 
    return false; }

	//===================================================== Reading every line and determining size
	int x = 0;
	int y = 0;
	int xmax = 0;
	int ymax = 0;
	while (feof(file)==0)
	{
		fgets(fileLine, 19999, file);    //Reads a single line
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
			if (x > xmax) xmax=x;
			if (y > ymax) ymax=y;
		}
	}
	rewind(file);

	//===================================================== Creates the matrix from the file
	this->SetSizeZero(ymax, xmax);
	x = 0;
	y = 0;
	while (feof(file) == 0)
	{
		fgets(fileLine, 19999, file);    //Reads a single line
		if (feof(file) == 0)
		{
			word = strtok(fileLine, " ");
			if (word != NULL)  {x = 1; y++;}
			while ((word != NULL) && (word[0] != char(10)))
			{
				sscanf(word, "%lf\n", &rawValue1);
				this->sij(y, x,rawValue1);
				x++;
				word = strtok(NULL, " ");
			}
			x--;
		}
	}
	fclose(file);

	return true;
}


//############################################################################# Create from lexical vector
/**This routine initialized the matrix from a lexicographically
sorted vector. The vector must be of size \f$ n^2 \f$ which
will result in the creation of a \f$ n\times n \f$ matrix.*/
void    CHI_MATH_MATRIX::CreateFromLexicalVector(CHI_MATH_MATRIX* x, long int size)
{
	//===================================================== Initializing entries
	if (!this->initialized)
	{
		this->SetSizeZero(size, size);
	}
	else
	{
		this->Resize(size, size);
	}
	
	//===================================================== Copying the values
	long int n = size;
	for (long int i = 1; i <= n; i++)
	{
		for (long int j = 1; j <= n; j++)
		{
			long int m = (i - 1)*n + j;
			this->sij(i, j,x->ij(m));
		}
	}



	return;
}


//############################################################################# Accessing a value
/** The entries of a matrix are private and therefore,
this routine is used to gain access to the values. Will return
a reference to \f$ a_{1,1} \f$.*/
double CHI_MATH_MATRIX::ij(long int row, long int col)
{
	if ((row>this->rowCount) || (col>this->colCount))
	{
		//MessageBox(NULL, "Invalid indexer", "Error", MB_ICONERROR);
		return this->entries[1][1];
	}
	return this->entries[row][col];
}

//############################################################################# Accessing a value
/** The entries of a matrix are private and therefore,
this routine is used to gain access to the values. Will return
a reference to \f$ a_{1,1} \f$.*/
void CHI_MATH_MATRIX::sij(long int row, long int col,const double value)
{
	if ((row>this->rowCount) || (col>this->colCount))
	{
		//MessageBox(NULL, "Invalid indexer", "Error", MB_ICONERROR);
		return;
	}
    this->entries[row][col]=value;
}
