#include "chi_math_linear_interpolation.h"
#include "chi_math_linear_indexer.h"



//############################################################################# Constructor
/** Constructor for indexer.*/
CHI_MATH_LINEAR_INDEXER::CHI_MATH_LINEAR_INDEXER()
{
	this->nextLevelDefined = false;
	minIndex=0;
	maxIndex=0;
	minValue=0.0;
	maxValue=0.0;
}

//############################################################################# Create recursion
/** Create a recursion of an indexer object.*/
void CHI_MATH_LINEAR_INDEXER::CreateRecursion(CHI_MATH_MATRIX* matrix, int levelAmount, int intervalsPerLevel)
{
	long int indexIntervals = (long int)(floor((this->maxIndex - this->minIndex) / intervalsPerLevel));

	long int min_counter = this->minIndex;
	long int max_counter = this->minIndex + indexIntervals - 1;
	for (int k = 1; k <= intervalsPerLevel; k++)
	{
		CHI_MATH_LINEAR_INDEXER* newIndexer = new CHI_MATH_LINEAR_INDEXER;
		newIndexer->minIndex = min_counter;
		newIndexer->maxIndex = max_counter;
		newIndexer->minValue = matrix->ij(min_counter, 1);
		newIndexer->maxValue = matrix->ij(max_counter, 1);
		newIndexer->level = this->level + 1;

		if (k == intervalsPerLevel)
		{
			newIndexer->maxIndex = this->maxIndex;
			newIndexer->maxValue = this->maxValue;
		}

		if (levelAmount > newIndexer->level)
		{
			newIndexer->nextLevelDefined = true;
			newIndexer->CreateRecursion(matrix, levelAmount, intervalsPerLevel);
		}

		this->indexers.PushItem(newIndexer);

		min_counter += indexIntervals;
		max_counter += indexIntervals;
		/*if (max_counter > this->maxIndex)
		{
			max_counter = this->maxIndex;
		}*/
	}

}


//############################################################################# Find indexes
/** This routine determines the internal indexes.*/
void     CHI_MATH_LINEAR_INDEXER::FindIndexRanges(double lookUpValue, long int& kmin, long int& kmax,CHI_MATH_MATRIX* matrix)
{
	minIndex=0;
	maxIndex=0;
	minValue=0.0;
	maxValue=0.0;
	if (!nextLevelDefined)
	{
		kmin = minIndex;
		kmax =   maxIndex;
		return;
	}

	CHI_MATH_LINEAR_INDEXER* indexer;
	CHI_MATH_LINEAR_INDEXER* lowrIndexer;
	CHI_MATH_LINEAR_INDEXER* highIndexer;
	//long int kkmin = 1;
	//long int kkmax = 1;
	//int indexColumn = 1;
	//int lookUpColumn = 2;
	//===================================================== Ruling out lowest value
	lowrIndexer = indexers.GetItem(0);
	if (lookUpValue <= lowrIndexer->minValue)
	{
		kmin = lowrIndexer->minIndex;
		kmax = kmin;
		return;
	}

	//===================================================== Ruling out the highest value
	highIndexer = indexers.GetItem(indexers.itemCount - 1);
	if (lookUpValue >= highIndexer->maxValue)
	{
		kmax = highIndexer->maxIndex;
		kmin = kmax;
		return;
	}

	//===================================================== Finding the indexer
	for (int k = 0; k < indexers.itemCount; k++)
	{
		indexer = indexers.GetItem(k);

		if ((lookUpValue >= indexer->minValue) && (lookUpValue <= indexer->maxValue))
		{
			break;
		}

	}

	//===================================================== Assigning scan indexes
	indexer->FindIndexRanges(lookUpValue, kmin, kmax,matrix);
}
