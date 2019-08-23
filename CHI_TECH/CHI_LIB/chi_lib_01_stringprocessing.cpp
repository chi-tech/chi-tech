#include "chi_lib.h"

#include <cstdio>
#include <cstdlib>

//#############################################################################
//############################### String comparison ###########################
/**Checks two null-terminated strings for equivalence.

\author Nakter*/
bool CHI_LIB::StringCompare(const char *input, const char *comparison)
{
	bool returnvalue=true;
	int		k;

	k=0;

	while (int(comparison[k]) != 0)
	{
		if (!(input[k] == comparison[k]))
		{
			returnvalue = false;
		}
		k+=1;
	}
	if (comparison[k] != input[k]) returnvalue = false;

	return returnvalue;
}

//#############################################################################
//############################### String copy #################################
/**Copies one string into another.

\param destination		Pointer to the destination.
\param source			Pointer to the source to be copied into the source.

\return					The location in the string of the null terminator.

\warning				This routine can produce heap corruption if source > destination.

\author Nakter*/
int CHI_LIB::StringCopy(char *destination, const char *source)
{
	int k=0;
	while (source[k] != '\0')
	{
		destination[k] = source[k];
		k++;
	}
	destination[k] = '\0';

	return k;
}

//#############################################################################
//############################ String concatenation ###########################
/**This function concatenates two strings to a newly creating string. The two
parameters are concatenated by removing Source1's trailing terminator and then
adding both these strings to a heap-allocated string.

\warning This function requires garbage collection since it creates heap memory
         without deleting it. If used, use with a pointer and delete after use.

\param Source1		Pointer to the first string to be concatenated.
\param Source2		Pointer to the second string to be concatenated.
\param sizeLimit	Limit the size of the newly created string;

\author Nakter*/
char* CHI_LIB::StringCat(char *Source1, char *Source2, int sizeLimit,bool deleteSource1)
{
	int		stringLength1 = CHI_LIB::GetStringLength(Source1);
	int		stringLength2 = CHI_LIB::GetStringLength(Source2);

	if (sizeLimit == -1)
	{
		sizeLimit = stringLength1+stringLength2+2;
	}

	char*	NewString = new char[sizeLimit];
	int		k,m,n;

	//===================================================== Copying Source1
	k=0;
	m=0;
	while ((Source1[m] != '\0') && (m<=stringLength1) && (m<sizeLimit) && (k<sizeLimit))
	{
		NewString[k] = Source1[m];
		k+=1;
		m+=1;
	}

	//===================================================== Copying Source2
	n=0;
	while ((Source2[n] != '\0') && (n<=stringLength2) && (n<sizeLimit) && (k<sizeLimit))
	{
		NewString[k] = Source2[n];
		k+=1;
		n+=1;
	}

	NewString[k] = '\0';

	//===================================================== Deleting the initial source if required
	if (deleteSource1) { delete [] Source1; }

	return NewString;

}

/**This routine receives the location of the Source1 null terminator and its capacity. If Source1
plus Source2 is greater than the capacity; the capacity is doubled and the string is copied into 
a new memory space.*/
char* CHI_LIB::StringCatEx2(char **Source1, char *Source2, int* currentSize, int* capacity)
{

	int		stringLength1 = *currentSize;
	int		stringLength2 = CHI_LIB::GetStringLength(Source2);

	//int sizeLimit = stringLength1+stringLength2+2;
	*currentSize = stringLength1+stringLength2;

	//===================================================== Checking if new object is required
	char*	NewString;
	if ((stringLength1+stringLength2+2) > *capacity) 
	{
		NewString = new char[*capacity*2]; 
		*capacity*=2;
		int k=0;
		for (k=0; k<=stringLength1; k++)
		{
			char* tempSource=*Source1;
			NewString[k] = tempSource[k];
		}
		NewString[k+1]='\0';
		delete [] *Source1;
	}
	//===================================================== If current capacity is sufficient
	else 
	{
		NewString = *Source1;
	}
	
	//===================================================== Adding Source2 to Source1
	int m=stringLength1;
	for (int k=0; k<=stringLength2; k++)
	{
		NewString[m] = Source2[k];
		m++;
	}
	//NewString[m]='\0';			//Adding a null terminator for good measure.

	////===================================================== Deleting the initial source if required
	//if (deleteSource1) { delete [] Source1; }
	*Source1 = NewString;
	return NewString;

}

//#############################################################################
//############################ String concatenation ###########################
/**This routine concatenates the second string to the first. This routine is a 
more memory-leak safe version of SmiLib::StringCat since it attaches the 
second string, Source2, to the first string, Source1.

\param Source1		Pointer to the first string.
\param Source2		Pointer to the second string.
\param sizeLimit	Buffer overrun size of the first string.

\author Nakter*/
void CHI_LIB::StringCatEx(char *Source1, char *Source2, int sizeLimit)
{
	int lineLength1 = CHI_LIB::GetStringLength(Source1);
	int lineLength2 = CHI_LIB::GetStringLength(Source2);

	int k=lineLength1;
	int m=0;
	while ((Source2[m] != '\0')&& (m<lineLength2)/* && (m<sizeLimit) && (k<(sizeLimit-1))*/)
	{
		Source1[k] = Source2[m];
		k+=1;
		m+=1;
	}
	Source1[k] = '\0';
}



//#############################################################################
//###################### GET AMOUNT OF LINES ##################################
/**Returns the amount of lines in a string. This function determines how many
lines are present in a string by examining the amount of char(10)s and 
char(13)s.

\param Pointer to a input string.

\return Amount of char(10)s plus char(13)s plus 1.

\author Nakter*/
int CHI_LIB::GetAmountOfLines(char *inputString)
{
	int linecount,k;										//Counter for the amount of lines
	
	//===================================================== Initiating run through all input characters
	k=0;
	linecount = 0;
	while (inputString[k] != '\0')
	{
		//================================================= Filtering end of a line
		if (k>0)
		{
			if ((inputString[k]==char(10)) || (inputString[k]==char(13)))
			{
				linecount+=1;
			}
		}
		k+=1;
	}
	linecount+=1;
	//=====

	//===================================================== Returning the word count
	return linecount;
}

//#############################################################################
//########################### GET STRING LENGTH ###############################
/**Gets location of first null terminator. This function is an indication of a 
strings size.

\param inputString	Pointer to the input string.
\param sizeLimit	Buffer over run position of the input string.

\author Nakter*/
int CHI_LIB::GetStringLength(char* inputString, int sizeLimit)
{
	int k=0;
	if (sizeLimit >= 0){	while ((inputString[k] != '\0') && (k < sizeLimit))	{k++;}}
	else			   {	while (inputString[k] != '\0')	{k++;}}
	return k;
}

//#############################################################################
//################################ PROCESS WORDS ##############################
/**Extracts a collection of words from an input string.

\author JIC Vermaak*/
int CHI_LIB::GetWords(char *InputString,CST_WORDCOLLECTION *WordCollection,int lineNumber)
{
	int k,m,n;																	//General use counters
	int linecount;																//Counter for the amount of lines
	int wordcount=0;															//Counter for the number of words
	int wordstart,wordend;														//Counters for word extraction
	char StdDelimiter=' ';														//Setting the standard delimiter
	
	//===================================================== Purging the word collection structure
	for (int a=0;a<100;a++)
	{
		CHI_LIB::StringCopy(WordCollection->Words[a]," ");
	}

	//===================================================== Initiating run through all input characters
	k=0;
	wordstart=0;
	wordend = -1;
	linecount = 0;
	while (InputString[k] != '\0')
	{
		//================================================= Filtering end of a line
		if (k>0)
		{
		if ((InputString[k]==char(10)) || (InputString[k]==char(13)))
		{
			linecount+=1;
			wordstart=k;
			wordend = k-1;
			if (linecount > lineNumber) break;
		}
		}

		//================================================= Filtering start of a word
		if (((InputString[k]==StdDelimiter) || (InputString[k]==char(9))) && ((InputString[k+1]!=StdDelimiter) && (InputString[k+1]!=char(9))))
		{
			wordstart = k+1;
		}

		//================================================= Filtering end of a word
		if (((InputString[k]!=StdDelimiter) && (InputString[k]!=char(9))) && ( (InputString[k+1]==StdDelimiter) || (InputString[k+1]==char(9)) || (InputString[k+1]==char(10)) || (InputString[k+1]==char(13)) || (InputString[k+1]=='\0') ))
		{
			wordend = k;
		}

		//================================================= Extracting a word
		if ((wordend >= wordstart) && (linecount == lineNumber))
		{
			wordcount += 1;
			//WordCollection->Words[wordcount] = new char[];
			n=-1;
			
			//============================================= Copying the word and adding a null terminator
			for (m=wordstart;m<=wordend;m++)
			{
				n+=1;
				WordCollection->Words[wordcount][n] = InputString[m];
			}
			WordCollection->Words[wordcount][n+1] = '\0';

			wordstart = k+1;
			wordend = k-1;
		}
		k+=1;
		
	}
	//=====

	//===================================================== Returning the word count
	if (wordcount >0)
	{
		return wordcount;
	}
	else return 0;
}




//#############################################################################
//###################### CONVERT NUM TO STRING ################################
/***/
char* CHI_LIB::StringFNum(int number)
{
	char* output = new char[20];

	//sprintf_s(output,20,"%d",number);
	sprintf(output,"%d",number);

	return output;
}

char* CHI_LIB::StringFNum(float number)
{
	char* output = new char[20];

	//sprintf_s(output,20,"%f",number);
	sprintf(output,"%f",number);

	return output;
}

char* CHI_LIB::StringFNum(float number,char* format)
{
	char* output = new char[20];

	sprintf(output, format, number);
	//sprintf_s(output, 20, format, number);

	return output;
}


//############################################################################# StringType
/**   Checks string to see what data types are stored and returns a number.
\author GMO*/
int CHI_LIB::StringType(char* inputString)
{
	int typeword = 0;
	char* endos;
	strtod(inputString, &endos);
	//int integer;

	if(!CHI_LIB::StringCompare("\0",endos))
	{
		//Return value for a string
		typeword = 1;
	}
	else
	{
		//integer = strtol(inputString,&endos,10);

		if(!CHI_LIB::StringCompare(endos,"\0"))
		{
			//Return value for a double
			typeword = 3;
		}
		else
		{
			//Return value for a integer
			typeword = 2;
		}
	}

	return typeword;
}







//############################################################################# String Insert
/**Inserts a string at the given position.
\author JIC*/
int CHI_LIB::StringInsert(char* original,int position,char* input)
{
	char temp[200];
	char temp2[200];

	for (int k=0;k<position;k++)
	{
		temp[k]=original[k];
	}
	temp[position]='\0';

	int length=CHI_LIB::GetStringLength(original);
	for (int k=position;k<length;k++)
	{
		temp2[k-position]=original[k];
	}
	temp2[length-position]='\0';

	original[0]='\0';
	CHI_LIB::StringCatEx(original,temp);
	CHI_LIB::StringCatEx(original,input);
	CHI_LIB::StringCatEx(original,temp2);

	return 0;
}





//############################################################################# Creates a new string that can be added
/**Creates a new string that can be added to a stack.
\author JIC*/
char* CHI_LIB::MakeString(char* stringContents)
{
	char* newString=new char[200];
	newString[0]='\0';

	CHI_LIB::StringCatEx(newString,stringContents);

	return newString;
}
