#ifndef CHI_TOOL_CLASS_VECTOR_H
#define CHI_TOOL_CLASS_VECTOR_H

#include <stdio.h>

//######################################################### CLASS DEFINITION
/**A class for the custom implementation of a vector data structure. This
template class can contain multiple polong inters to any structure, class or 
variable, of any type. Items are added and removed in a strictly Object-
Orientated manner to avoid buffer overruns. The vector can also exhibit a
stack type behaviour with the pushing and popping of items.

\param VectorType The type specification for the item to be stored.

\author Nakter*/
template<class VectorType>
class CHI_VECTOR
{
	typedef VectorType* VectorTypePtr;				///< Polong inter typedefinition for template class parameter.
	//===================================================== Public member variables
public:
	long int		capacity;											///< Current vector capacity.
	long int		expansionFactor;							///< Factor by which capacity is increased when reached.
	long int		itemCount;										///< Amount of items in the vector.
	long int		stackCount;										///< Amount of array positions from top to bottom, filled or not.
	bool        threadProtectionEnabled;      ///< Setting to allow thread protection (Default: True)
private:
	VectorTypePtr	*Item;											///< Dynamically allocated array of items.

public:
	//===================================================== Public member functions
	//===== Constructor
					CHI_VECTOR();

	long int			AddItem(VectorTypePtr inputItem);
	long int			PushItem(VectorTypePtr inputItem);
	void					SetItem(long int index, VectorTypePtr inputItem);
	long int			CreateAddItem();
	long int			CreatePushItem();
	VectorTypePtr GetItem(long int index);
	VectorTypePtr GetFirstNonNULLItem();
	VectorTypePtr	GetPopItem();						
	VectorTypePtr	PopItem();
	VectorTypePtr	PullItem(long int index);
	VectorTypePtr SwapItem(long int index, VectorTypePtr inputItem);
	void          InsertItem(long int index,VectorTypePtr inputItem);
	void					KickItem(long int index);
	void					ClearItem(long int index);
	void					ClearVector();
	void					EmptyVector();
};

//#############################################################################
//######################### DEFAULT CONSTRUCTOR ###############################
/**Default constructor for the CHI_VECTOR class.

\version NVC
\author Nakter*/
template<class VectorType>
CHI_VECTOR<VectorType>::CHI_VECTOR()
{
	//======================= Initializing member variables
	capacity = 100;
	expansionFactor = 2;
	itemCount = 0;
	stackCount = 0;
	threadProtectionEnabled = true;

	//======================= Creating array of polong inters
	Item = new VectorTypePtr[capacity];

	
	for (long int k=0;k<capacity;k++) {Item[k]=NULL;}
}

//#############################################################################
//############################# ADD ITEM ######################################
/**Method for adding an item to the vector. This method adds an item of 
VectorType to the back of the vector. The argument given is a polong inter to an 
ALREADY created object.

\param inputItem A polong inter to the item to be added to the vector.

\version 1.0 2011-07-07
\author Nakter*/
template<class VectorType>
long int CHI_VECTOR<VectorType>::AddItem(VectorTypePtr inputItem)
{
	
	long int emptyIndex = -1;

	if (inputItem == NULL)
	{
		itemCount -= 1;
		return emptyIndex;
	}

	//===================================================== Find possible long interstitial space and return if found
	for (long int k=0; k<stackCount; k++)
	{
		if (Item[k] == NULL)
		{
			Item[k]		= inputItem;
			emptyIndex	= k;
			return emptyIndex;
		}
	}

	//===================================================== No long intersticial space
	itemCount += 1;
	stackCount += 1;

	//===================================================== Adjusting capacity
	if (stackCount > capacity)
	{
		//=================== Making temporary copy
		VectorTypePtr *tempVector = Item;

		//=================== Rebuilding new array
		Item = new VectorTypePtr[capacity*expansionFactor];										//New larger array
		for (long int k=0;k<capacity;k++)							{ Item[k] = tempVector[k]; }	//Copying existing items
		for (long int k=capacity;k<capacity*expansionFactor;k++)		{ Item[k] = NULL; }				//Clear new spaces

		//=================== Deleting temp vector
		delete [] tempVector;

		//=================== Updating the capacity
		capacity = capacity*expansionFactor;
	}

	//===================================================== Assigning the polong inter
	Item[stackCount-1] = inputItem;

	return (stackCount-1);
}

//#############################################################################
//############################# PUSH ITEM #####################################
/**Method for adding an item to the top of the vector. This method adds an item of 
VectorType to the top of the vector. The argument given is a polong inter to an 
ALREADY created object.

\param inputItem A polong inter to the item to be added to the vector.

\version 1.0 2011-07-07
\author Nakter*/
template<class VectorType>
long int CHI_VECTOR<VectorType>::PushItem(VectorTypePtr inputItem)
{
	itemCount += 1;
	stackCount += 1;

	//===================================================== Adjusting capacity
	if (stackCount > capacity)
	{
		//=================== Making temporary copy
		VectorTypePtr *tempVector = Item;

		//=================== Rebuilding new array
		Item = new VectorTypePtr[capacity*expansionFactor];										//New larger array
		for (long int k=0;k<capacity;k++)							{ Item[k] = tempVector[k]; }	//Copying existing items
		for (long int k=capacity;k<capacity*expansionFactor;k++)		{ Item[k] = NULL; }				//Clear new spaces

		//=================== Deleting temp vector
		delete [] tempVector;

		//=================== Updating the capacity
		capacity = capacity*expansionFactor;
	}

	//===================================================== Assigning the polong inter
	Item[stackCount-1] = inputItem;

	return (stackCount-1);
}




template<class VectorType>
void CHI_VECTOR<VectorType>::SetItem(long int index, VectorTypePtr inputItem)
{
	//===================================================== Check query within stack count bounds
	if ((index > this->stackCount) || (index < 0)) {return;}

	//===================================================== Inserting the item
	this->Item[index]=inputItem;

}

//#############################################################################
//####################### CREATE ADD ITEM #####################################
/**Method for adding an item to the vector by creating it first. This method 
adds an item of VectorType to the back of the vector.

\param inputItem A polong inter to the item to be added to the vector.

\version 1.0 2011-07-07
\author Nakter*/
template<class VectorType>
long int CHI_VECTOR<VectorType>::CreateAddItem()
{
	VectorTypePtr newItem = new VectorType;
	itemCount += 1;
	long int emptyIndex = -1;

	//===================================================== Find possible long interstitial space and return if found
	for (long int k=0; k<stackCount; k++)
	{
		if (Item[k] == NULL)
		{
			Item[k]		= newItem;
			emptyIndex	= k;
			return emptyIndex;
		}
	}

	//===================================================== No long intersticial space
	stackCount += 1;

	//===================================================== Adjusting capacity
	if (stackCount > capacity)
	{
		//=================== Making temporary copy
		VectorTypePtr *tempVector = Item;

		//=================== Rebuilding new array
		Item = new VectorTypePtr[capacity*expansionFactor];										//New larger array
		for (long int k=0;k<capacity;k++)							{ Item[k] = tempVector[k]; }	//Copying existing items
		for (long int k=capacity;k<capacity*expansionFactor;k++)		{ Item[k] = NULL; }				//Clear new spaces

		//=================== Deleting temp vector
		delete [] tempVector;

		//=================== Updating the capacity
		capacity = capacity*expansionFactor;
	}

	//===================================================== Assigning the polong inter
	Item[stackCount-1] = newItem;

	return (stackCount-1);
}

//#############################################################################
//####################### CREATE PUSH ITEM ####################################
/**Method for adding an item to the top of the vector by creating it first. 
This method adds an item of VectorType to the top of the vector. The argument 
given is a polong inter to an ALREADY created object.

\param inputItem A polong inter to the item to be added to the vector.

\version 1.0 2011-07-07
\author Nakter*/
template<class VectorType>
long int CHI_VECTOR<VectorType>::CreatePushItem()
{
	VectorTypePtr newItem = new VectorType;
	itemCount += 1;
	stackCount += 1;

	//===================================================== Adjusting capacity
	if (stackCount > capacity)
	{
		//=================== Making temporary copy
		VectorTypePtr *tempVector = Item;

		//=================== Rebuilding new array
		Item = new VectorTypePtr[capacity*expansionFactor];										//New larger array
		for (long int k=0;k<capacity;k++)							{ Item[k] = tempVector[k]; }	//Copying existing items
		for (long int k=capacity;k<capacity*expansionFactor;k++)		{ Item[k] = NULL; }				//Clear new spaces

		//=================== Deleting temp vector
		delete [] tempVector;

		//=================== Updating the capacity
		capacity = capacity*expansionFactor;
	}

	//===================================================== Assigning the polong inter
	Item[stackCount-1] = newItem;

	return (stackCount-1);
}

//#############################################################################
//########################## GET AN ITEM ###################################
/**Method for getting the specified item. This routine returns the required
polong inter but only does so if the index is valid and the item exists. The design
choice was to return NULL if an index is required beyond the stackCount value.
This can be re-evaluated if better feedback is required.

\param index	Location in the stack of the vector where the polong inter is.

\author Nakter*/
template<class VectorType>
VectorType* CHI_VECTOR<VectorType>::GetItem(long int index)
{
	//===================================================== Check query within stack count bounds
	if ((index > this->stackCount) || (index < 0)) {return NULL;}

	//===================================================== If item is non-null return polong inter
	if (Item[index] != NULL)
	{
		return Item[index];
	}

	//===================================================== If item is null return null
	else
	{
		return NULL;
	}
}

//#############################################################################
//########################## GET AN ITEM ###################################
/**Method for getting the specified item. This routine returns the required
polong inter but only does so if the index is valid and the item exists. The design
choice was to return NULL if an index is required beyond the stackCount value.
This can be re-evaluated if better feedback is required.

\param index	Location in the stack of the vector where the polong inter is.

\author Nakter*/
template<class VectorType>
VectorType* CHI_VECTOR<VectorType>::GetFirstNonNULLItem()
{
	for (long int index=0;index<stackCount;index++)
	{
		//===================================================== Check query within stack count bounds
		if ((index > this->stackCount) || (index < 0)) {return NULL;}
		
		//===================================================== If item is non-null return polong inter
		if (Item[index] != NULL)
		{
			return Item[index];
		}
			
			//===================================================== If item is null return null
		else
		{
			return NULL;
		}
	}
	return NULL;
}

//#############################################################################
//########################## GET POP AN ITEM ###################################
/**Method for getting the top item without removing it. This routine returns the 
top polong inter but only does so if the item exists.

\param index	Location in the stack of the vector where the polong inter is.

\author Nakter*/
template<class VectorType>
VectorType* CHI_VECTOR<VectorType>::GetPopItem()
{

	//===================================================== Check query within stack count bounds
	if (stackCount <= 0) { return NULL; }

	//===================================================== If item is non-null, return polong inter
	if (Item[stackCount-1] != NULL)
	{
		return Item[stackCount-1];
	}

	//===================================================== If item is null, return null
	else
	{
		return NULL;
	}
}

//#############################################################################
//############################# POP AN ITEM ###################################
/**Method for getting the top item by removing it. This routine returns the 
top polong inter but only does so if the item exists.

\param index	Location in the stack of the vector where the polong inter is.

\author Nakter*/
template<class VectorType>
VectorType* CHI_VECTOR<VectorType>::PopItem()
{
	VectorTypePtr returnValue;

	//===================================================== Check query within stack count bounds
	if (stackCount <= 0) {return NULL; }

	//===================================================== If item is non-null, return polong inter
	if (Item[stackCount-1] != NULL)
	{
		returnValue =  Item[stackCount-1];

		//============================= Clearing item
		Item[stackCount-1] = NULL;
		this->itemCount	-= 1;
		this->stackCount -= 1;
	}

	//===================================================== If item is null, return null
	else
	{
		return NULL;
	}

	return returnValue;
}

//#############################################################################
//########################## PULL AN ITEM #####################################
/**Removes a given polong inter, without deleting its value and drops the vector to 
assimilate the gap.

\author JIC Vermaak*/
template<class VectorType>
VectorType* CHI_VECTOR<VectorType>::PullItem(long int index)
{
	VectorTypePtr returnValue;
	//===================================================== Check query within stack count bounds
	if (index > this->stackCount) {return NULL;}

	//===================================================== Get the item
	if (this->Item[index] != NULL)
	{
		returnValue = this->Item[index];
	}
	else
	{
		return NULL;
	}
	
	//===================================================== Assimilate the gap
	for (long int k=index; k<(this->itemCount -1);k++)
	{
		this->Item[k] = this->Item[k+1];
	}
	Item[itemCount-1] = NULL;
	itemCount -= 1;
	stackCount-=1;

	return returnValue;
}

//#############################################################################
//####################### Hotswap Item ####################################
/** Method for hot swapping a given index from the source vector to the target vector.
 * This method will only swap the index if the two vectors indexes match.

\param index Index to be swapped from source to target
\param inputItem Target Vector

\version 1.0 2017/09/06
\author Gmo*/
template<class VectorType>
VectorType* CHI_VECTOR<VectorType>::SwapItem(long int index, VectorTypePtr inputItem)
{
	VectorTypePtr tempVector = NULL;

	//===================================================== Check query within stack count bounds
	if((index > this->stackCount) || (index < 0)){return NULL;};

	tempVector = Item[index];
	Item[index] = inputItem;

	return tempVector;
}

//#############################################################################
//######################### INSERT AN ITEM ####################################
/**Inserts the item at the given index and shifts the items up.

\author JIC Vermaak*/
template<class VectorType>
void CHI_VECTOR<VectorType>::InsertItem(long int index,VectorTypePtr inputItem)
{
	//===================================================== Check query within stack count bounds
	if (index > this->stackCount) { return; }

	//===================================================== Pushing the last item on top of the stack
	long int last_size = this->itemCount;
	this->PushItem(Item[last_size-1]);


	//===================================================== Push all items up
	for (long int k=last_size-2; k>=index;k--)
	{
		this->Item[k+1] = this->Item[k];
	}

	//===================================================== Insert new item
	Item[index] = inputItem;
}

//#############################################################################
//########################## KICK AN ITEM #####################################
/**Removes a given polong inter, deletes its value and drops the vector to assimilate
the gap.

\author JIC Vermaak*/
template<class VectorType>
void CHI_VECTOR<VectorType>::KickItem(long int index)
{
	//===================================================== Check query within stack count bounds
	if (index > this->stackCount) { return; }

	//===================================================== Delete the item
	if (this->Item[index] != NULL)
	{
		delete [] this->Item[index];
	}
	else
	{
		return;
	}
	
	//===================================================== Assimilate the gap
	for (long int k=index; k<(this->itemCount -1);k++)
	{
		this->Item[k] = this->Item[k+1];
	}
	Item[itemCount-1] = NULL;
	itemCount -= 1;
	stackCount -= 1;
}

//#############################################################################
//########################## CLEAR AN ITEM ###################################
/**Method for clearing a polong inter from the vector. This routine 
removes a polong inter from the vector stack at the indicated index by setting it to
NULL. It does not assimilate the empty spot.

\param index	Location in the stack of the vector where the polong inter is.

\author Nakter*/
template<class VectorType>
void CHI_VECTOR<VectorType>::ClearItem(long int index)
{

	//===================================================== Check query within stack count bounds
	if (index > this->stackCount) { return; }

	//=======================================
	Item[index] = NULL;
}

//#############################################################################
//############################ CLEAR THE VECTOR ###############################
/**Method for emptying the entire vector stack.

\author Nakter*/
template<class VectorType>
void CHI_VECTOR<VectorType>::ClearVector()
{
	for (long int k=(itemCount-1); k>=0; k--)
	{
		if (Item[k] != NULL)
		{
		delete [] Item[k];
		Item[k] = NULL;
		}
	}

	itemCount = 0;
	stackCount = 0;
}

//#############################################################################
//############################ EMPTY THE VECTOR ###############################
/**Method for emptying the entire vector stack but not deleting anything.

\author Nakter*/
template<class VectorType>
void CHI_VECTOR<VectorType>::EmptyVector()
{
	for (long int k=(itemCount-1); k>=0; k--)
	{
		if (Item[k] != NULL)
		{

		Item[k] = NULL;
		}
	}

	itemCount = 0;
	stackCount = 0;
}

#endif
