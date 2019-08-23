#include"chi_math_linear_interpolation.h"




//############################################################################# Constructor
/** Default constructor*/
CHI_MATH_LINEAR_INTERPOLATION::CHI_MATH_LINEAR_INTERPOLATION()
{
	this->intervalsPerLevel = 10;
	this->numberOfLevels = 5;
}

//############################################################################# Linear interpolation Optimized
/** Makes an optimized linear interpolation given a lookup value and a lookup table.*/
double CHI_MATH_LINEAR_INTERPOLATION::LinearInterpolate(double lookUpValue)
{
	CHI_MATH_LINEAR_INDEXER* indexer;
	CHI_MATH_LINEAR_INDEXER* lowrIndexer;
	CHI_MATH_LINEAR_INDEXER* highIndexer;
	long int kmin = 1;
	long int kmax = 1;
	int indexColumn = 1;
	int lookUpColumn = 2;
	//===================================================== Ruling out lowest value
	lowrIndexer = this->indexers.GetItem(0);
	if (lookUpValue <= lowrIndexer->minValue)
	{
		kmin = lowrIndexer->minIndex;
		return this->mainMatrix->ij(kmin, lookUpColumn);
	}

	//===================================================== Ruling out the highest value
	highIndexer = this->indexers.GetItem(this->indexers.itemCount-1);
	if (lookUpValue >= highIndexer->maxValue)
	{
		kmax = highIndexer->maxIndex;
		return this->mainMatrix->ij(kmax, lookUpColumn);
	}

	//===================================================== Finding the indexer
	indexer = this->indexers.GetItem(0);

	for (int k = 0; k < this->indexers.itemCount; k++)
	{
		indexer = this->indexers.GetItem(k);

		if ((lookUpValue >= indexer->minValue) && (lookUpValue <= indexer->maxValue))
		{
			break;
		}

	}

	//===================================================== Assigning scan indexes
	indexer->FindIndexRanges(lookUpValue, kmin, kmax,mainMatrix);


	//===================================================== If lookupValue <= smallest index
	if (lookUpValue <= this->mainMatrix->ij(kmin, indexColumn))
	{
		return this->mainMatrix->ij(kmin, lookUpColumn);
	}

	//===================================================== If lookupValue >= largest index
	if (lookUpValue >= this->mainMatrix->ij(kmax, indexColumn))
	{
		return this->mainMatrix->ij(kmax, lookUpColumn);
	}

	//===================================================== Handling linear interpolation
	for (long int k = kmin; k < kmax; k++)
	{
		if ((lookUpValue >= this->mainMatrix->ij(k, indexColumn)) && (lookUpValue < this->mainMatrix->ij(k + 1, indexColumn)))
		{

			double returnValue;
			returnValue = (this->mainMatrix->ij(k + 1, lookUpColumn) - this->mainMatrix->ij(k, lookUpColumn));
			returnValue = returnValue / (this->mainMatrix->ij(k + 1, indexColumn) - this->mainMatrix->ij(k, indexColumn));
			returnValue = returnValue * (lookUpValue - this->mainMatrix->ij(k, indexColumn));
			returnValue = returnValue + this->mainMatrix->ij(k, lookUpColumn);
			return returnValue;
		}
	}

	//===================================================== Problem in lookup Value

	return 0.0;
}

//############################################################################# Optimize index searches
/** This routine creates multiple levels of index searches.*/
void CHI_MATH_LINEAR_INTERPOLATION::Optimize(CHI_MATH_MATRIX* matrix)
{
	this->mainMatrix = matrix;
	long int indexIntervals = (long int)(floor(matrix->rowCount / this->intervalsPerLevel));

	long int min_counter=1;
	long int max_counter = indexIntervals;
	for (int k = 1; k <= this->intervalsPerLevel; k++)
	{
		CHI_MATH_LINEAR_INDEXER* newIndexer = new CHI_MATH_LINEAR_INDEXER;
		newIndexer->minIndex = min_counter;
		newIndexer->maxIndex = max_counter;
		newIndexer->minValue = matrix->ij(min_counter, 1);
		newIndexer->maxValue = matrix->ij(max_counter, 1);
		newIndexer->level = 1;

		if (k == this->intervalsPerLevel)
		{
			newIndexer->maxIndex = matrix->rowCount;
			newIndexer->maxValue = matrix->ij(matrix->rowCount, 1);
		}

		if (this->numberOfLevels > 1)
		{
			newIndexer->nextLevelDefined = true;
			newIndexer->CreateRecursion(matrix, this->numberOfLevels, this->intervalsPerLevel);
		}

		this->indexers.PushItem(newIndexer);

		min_counter += indexIntervals;
		max_counter += indexIntervals;
		/*if (max_counter > matrix->rowCount)
		{
			max_counter = matrix->rowCount;
		}*/
	}
}



