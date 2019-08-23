#include"chi_math_spline.h"
#include"../chi_math.h"

extern CHI_MATH                chi_math_handler;


//#pragma warning(disable: 4996)
//############################################################################# Initialize form file
/** Initializes a spline from file.

The file format follows the following logic.\n
\code
type[int]
minBoundary[double] maxBoundary[double]
Knot0[double]       KnotValue0[double]
Knot1[double]       KnotValue1[double]
Knot2[double]       KnotValue2[double]
Knot3[double]       KnotValue3[double]
etc.
\endcode
\return Returns false if failed.*/
bool CHI_MATH_SPLINE::InitializeFromFile(char* fileName)
{
	FILE* iFile;
	iFile = fopen(fileName,"r");
	int linesRead=0;

	double	value1=0.0;
	double	value2=0.0;
	int		splineType=0;
	double* pvalue1;
	double* pvalue2;

	if (iFile == NULL){ return false; }

	while (feof(iFile) == 0)
	{
		if (linesRead==0)
		{
			fscanf(iFile,"%d\n",&splineType);
			this->type=splineType;
		}
		fscanf(iFile,"%lf %lf\n",&value1,&value2);
		linesRead++;
		
		//================ Reading boundary conditions on top of file
		if (linesRead==1)
		{
			this->minBoundary=value1;
			this->maxBoundary=value2;
		}
		//================ Reading rest of values
		else
		{
			pvalue1=new double; *pvalue1=value1;
			pvalue2=new double; *pvalue2=value2;
			this->knotPoints.PushItem(pvalue1);
			this->knotValues.PushItem(pvalue2);
		}
	}

	fclose(iFile);

	//===================================================== Building the spline
	if (this->type==1){this->BuildNaturalSpline();}
	if (this->type==2){this->BuildLeftFixedSpline();}
	if (this->type==3){this->BuildRightFixedSpline();}
	if (this->type==4){this->BuildFixedSpline();}


	return true;
}

//############################################################################# Create from table
/**Initializes the spline given the table of knots and knot values.
\return Will always return true.*/
bool CHI_MATH_SPLINE::InitializeFromTable(CHI_MATH_MATRIX* table,int type,double minBound,double maxBound)
{
	//===================================================== Loading knots
	double* pvalue1;
	double* pvalue2;
	for (int i=1;i<=table->rowCount;i++)
	{
		pvalue1=new double; *pvalue1=table->ij(i,1);
		pvalue2=new double; *pvalue2=table->ij(i,2);
		this->knotPoints.PushItem(pvalue1);
		this->knotValues.PushItem(pvalue2);
	}

	//===================================================== Setting boundary conditions
	this->type=type;
	this->minBoundary=minBound;
	this->maxBoundary=maxBound;
	

	//===================================================== Building the spline
	if (this->type==1){this->BuildNaturalSpline();}
	if (this->type==2){this->BuildLeftFixedSpline();}
	if (this->type==3){this->BuildRightFixedSpline();}
	if (this->type==4){this->BuildFixedSpline();}


	return true;
}

//############################################################################# Build spline
/** Build the spline unknowns assuming free ends.*/
void CHI_MATH_SPLINE::BuildNaturalSpline()
{
	//===================================================== Initializing
	int n=this->knotPoints.itemCount-1;
	CHI_MATH_MATRIX A(n+1,n+1);
	CHI_MATH_MATRIX b(n+1,1);

	//===================================================== Filling the matrix
	double him1=0.0;
	double hi=0.0;
	double yim1=0.0;
	double yip1=0.0;
	double yi=0.0;
	for (int i=0;i<=(n);i++)
	{
		if(i<n) {hi  =*this->knotPoints.GetItem(i+1)-*this->knotPoints.GetItem(i);}
		if(i>0) {him1=*this->knotPoints.GetItem(i)-*this->knotPoints.GetItem(i-1);}

		if(i>0) {yim1=*this->knotValues.GetItem(i-1);}
		if(i<n) {yip1=*this->knotValues.GetItem(i+1);}
		yi  =*this->knotValues.GetItem(i);

		if(i==0)
		{
			A.sij(i+1,i+1,1.0);
			b.sij(i+1,1,0.0);
		}
		else if (i==n)
		{
			A.sij(i+1,i+1,1.0);
			b.sij(i+1,1,0.0);
		}
		else
		{
			A.sij(i+1,i  ,him1);

			A.sij(i+1,i+1  ,2.0*(hi+him1));

			A.sij(i+1,i+2,hi);

			b.sij(i+1,1,(6.0/hi)*(yip1-yi)-(6.0/him1)*(yi-yim1));
		}
	}

	//===================================================== Solve
	chi_math_handler.SolveUsingJacobiIteration(&A,this->z,&b);
}

//############################################################################# Build spline
/** Build the spline unknowns assuming a left fixed end.*/
void CHI_MATH_SPLINE::BuildLeftFixedSpline()
{
	//===================================================== Initializing
	int n=this->knotPoints.itemCount-1;
	CHI_MATH_MATRIX A(n+1,n+1);
	CHI_MATH_MATRIX b(n+1,1);

	//===================================================== Filling the matrix
	double him1=0.0;
	double hi=0.0;
	double yim1=0.0;
	double yip1=0.0;
	double yi=0.0;
	for (int i=0;i<=(n);i++)
	{
		if(i<n) {hi  =*this->knotPoints.GetItem(i+1)-*this->knotPoints.GetItem(i);}
		if(i>0) {him1=*this->knotPoints.GetItem(i)-*this->knotPoints.GetItem(i-1);}

		if(i>0) {yim1=*this->knotValues.GetItem(i-1);}
		if(i<n) {yip1=*this->knotValues.GetItem(i+1);}
		yi  =*this->knotValues.GetItem(i);

		if(i==0)
		{
			A.sij(i+1,i+1  ,2.0*hi);
			A.sij(i+1,i+2,hi);
			b.sij(i+1,1,(6.0/hi)*(yip1-yi)-6.0*this->minBoundary);
		}
		else if (i==n)
		{
			A.sij(i+1,i+1,1.0);
			b.sij(i+1,1,0.0);
		}
		else
		{
			A.sij(i+1,i,him1);

			A.sij(i+1,i+1  ,2.0*(hi+him1));

			A.sij(i+1,i+2,hi);

			b.sij(i+1,1,(6.0/hi)*(yip1-yi)-(6.0/him1)*(yi-yim1));
		}
	}

	//===================================================== Solve
	chi_math_handler.SolveUsingJacobiIteration(&A,this->z,&b);
}




//############################################################################# Build spline
/** Build the spline unknowns assuming a right fixed end.*/
void CHI_MATH_SPLINE::BuildRightFixedSpline()
{
	//===================================================== Initializing
	int n=this->knotPoints.itemCount-1;
	CHI_MATH_MATRIX A(n+1,n+1);
	CHI_MATH_MATRIX b(n+1,1);

	//===================================================== Filling the matrix
	double him1=0.0;
	double hi=0.0;
	double yim1=0.0;
	double yip1=0.0;
	double yi=0.0;
	for (int i=0;i<=(n);i++)
	{
		if(i<n) {hi  =*this->knotPoints.GetItem(i+1)-*this->knotPoints.GetItem(i);}
		if(i>0) {him1=*this->knotPoints.GetItem(i)-*this->knotPoints.GetItem(i-1);}

		if(i>0) {yim1=*this->knotValues.GetItem(i-1);}
		if(i<n) {yip1=*this->knotValues.GetItem(i+1);}
		yi  =*this->knotValues.GetItem(i);

		if(i==0)
		{
			A.sij(i+1,i+1,1.0);
			b.sij(i+1,1,0.0);
			
		}
		else if (i==n)
		{
			A.sij(i+1,i+1  ,2.0*him1);
			A.sij(i+1,i  ,him1);
			b.sij(i+1,1,(-1.0*6.0/him1)*(yi-yim1)+6.0*this->maxBoundary);
		}
		else
		{
			A.sij(i+1,i,him1);

			A.sij(i+1,i+1  ,2.0*(hi+him1));

			A.sij(i+1,i+2,hi);

			b.sij(i+1,1,(6.0/hi)*(yip1-yi)-(6.0/him1)*(yi-yim1));
		}
	}

	//===================================================== Solve
	chi_math_handler.SolveUsingJacobiIteration(&A,this->z,&b);
}



//############################################################################# Build spline
/** Build the spline unknowns assuming a fixed ends on both sides.*/
void CHI_MATH_SPLINE::BuildFixedSpline()
{
	//===================================================== Initializing
	int n=this->knotPoints.itemCount-1;
	CHI_MATH_MATRIX A(n+1,n+1);
	CHI_MATH_MATRIX b(n+1,1);

	//===================================================== Filling the matrix
	double him1=0.0;
	double hi=0.0;
	double yim1=0.0;
	double yip1=0.0;
	double yi=0.0;
	for (int i=0;i<=(n);i++)
	{
		if(i<n) {hi  =*this->knotPoints.GetItem(i+1)-*this->knotPoints.GetItem(i);}
		if(i>0) {him1=*this->knotPoints.GetItem(i)-*this->knotPoints.GetItem(i-1);}

		if(i>0) {yim1=*this->knotValues.GetItem(i-1);}
		if(i<n) {yip1=*this->knotValues.GetItem(i+1);}
		yi  =*this->knotValues.GetItem(i);

		if(i==0)
		{
			A.sij(i+1,i+1  ,2.0*hi);
			A.sij(i+1,i+2,hi);
			b.sij(i+1,1,(6.0/hi)*(yip1-yi)-6.0*this->minBoundary);
			
		}
		else if (i==n)
		{
			A.sij(i+1,i+1  ,2.0*him1);
			A.sij(i+1,i  ,him1);
			b.sij(i+1,1,(-1.0*6.0/him1)*(yi-yim1)+6.0*this->maxBoundary);
		}
		else
		{
			A.sij(i+1,i,him1);

			A.sij(i+1,i+1  ,2.0*(hi+him1));

			A.sij(i+1,i+2,hi);

			b.sij(i+1,1,(6.0/hi)*(yip1-yi)-(6.0/him1)*(yi-yim1));
		}
	}

	//===================================================== Solve
	chi_math_handler.SolveUsingJacobiIteration(&A,this->z,&b);
}








//############################################################################# Spline interpolation
/** Finds the spline interpolation value for x.

\todo Optimize for large splines (many knots).*/
double CHI_MATH_SPLINE::Value(double x)
{
	double hi=0.0;
	double yip1=0.0;    
	double yi=0.0;
	double zip1=0.0;	//
	double zi=0.0;      //
	double tip1=0.0;
	double ti=0.0;
	double n=this->knotPoints.itemCount-1;
	for (int i=0;i<=(n-1);i++)
	{
		if ((x>=*this->knotPoints.GetItem(i))&&(x<=*this->knotPoints.GetItem(i+1)))
		{
			zi=this->z->ij(i+1);
			zip1=this->z->ij(i+2);

			ti  =*this->knotPoints.GetItem(i);
			tip1=*this->knotPoints.GetItem(i+1);

			yi  =*this->knotValues.GetItem(i);
			yip1=*this->knotValues.GetItem(i+1);

			hi=tip1-ti;

			return (zi/6.0/hi)*pow(tip1-x,3.0) + (zip1/6.0/hi)*pow(x-ti,3.0) +
				   (yip1/hi-zip1*hi/6.0)*(x-ti) + (yi/hi - zi*hi/6.0)*(tip1-x);
		}
	}
	return *this->knotValues.GetItem(0);
}

//############################################################################# Spline interpolation
/** Finds the spline interpolation derivative value for x.

\todo Optimize for large splines (many knots).*/
double CHI_MATH_SPLINE::Derivative(double x)
{
	double hi=0.0;
	double yip1=0.0;    
	double yi=0.0;
	double zip1=0.0;	//
	double zi=0.0;      //
	double tip1=0.0;
	double ti=0.0;
	double n=this->knotPoints.itemCount-1;
	for (int i=0;i<=(n-1);i++)
	{
		if ((x>=*this->knotPoints.GetItem(i))&&(x<=*this->knotPoints.GetItem(i+1)))
		{
			zi=this->z->ij(i+1);
			zip1=this->z->ij(i+2);

			ti  =*this->knotPoints.GetItem(i);
			tip1=*this->knotPoints.GetItem(i+1);

			yi  =*this->knotValues.GetItem(i);
			yip1=*this->knotValues.GetItem(i+1);

			hi=tip1-ti;

			return (-1.0*zi/2.0/hi)*pow(tip1-x,2.0) + (zip1/2.0/hi)*pow(x-ti,2.0) +
				   (yip1/hi-zip1*hi/6.0) - (yi/hi - zi*hi/6.0);
		}
	}
	return *this->knotValues.GetItem(0);
}
