#ifndef CHI_LIB_H
#define CHI_LIB_H

#include "chi_lib_structs.h"

// ############################################################################# Namespace definition
/** Namespace
//\author G-mo */
namespace CHI_LIB
{
	//01 String processing
	bool					StringCompare(const char *input, const char *comparison);
	int						StringCopy(char* destination,const char* source);
	char*					StringCat(char *Source1, char *Source2,int sizeLimit=1000, bool deleteSource1=false);
	char*					StringCatEx2(char **Source1, char *Source2,int* currentSize=0,int* capacity=0);
	void	    	     	StringCatEx(char *Source1, char *Source2, int sizeLimit=1000);
	int						GetAmountOfLines(char* inputString);
	int						GetStringLength(char* inputString,int sizeLimit=-1);
	int						GetWords(char *InputString,CST_WORDCOLLECTION *WordCollection,int lineNumber=0);
	char*					StringFNum(int number);
	char*					StringFNum(float number);
	char*					StringFNum(float number,char* format);
	int                     StringType(char* inputString);
	int                     StringInsert(char* original,int position,char* input);
	char*                   MakeString(char* stringContents);

	//02 File control
	int 					TextFileWrite(char *filename, char *writeString, const char* option);
	char* 					TextFileRead(char *fileName);

	//03 Folder control
	bool 					CreateDirectoryS(char* sPath,bool verifyCompletion=true);
	bool 					DeleteDirectory(char* sPath);
	bool 					DeleteDirectoryContents(char* sPath);
	int 					GetListOfFiles(char* sPath,char* fileList[]);
};

#endif
