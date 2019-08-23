#include "chi_lib.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>


//############################################################################# WRITE A TEXT FILE
/**Writes the given string to a text file.

\param filename		Name of the file to be written.
\param writeString	Character string to write to the file.

return	Returns 1 if success and zero otherwise.

\author JIC Vermaak*/
int CHI_LIB::TextFileWrite(char *filename, char *writeString, const char* option)
{
	FILE		*filepointer;													//Pointer to the opened file
	//int			error;															//Error returned
	int			status = 0;														//Status of the opened file "1" indicates a successful write

	if (filename != NULL)														//If the file is not specified
	{
		//error = fopen_s(&filepointer,filename,option);								//Opening the file for overwriting
		filepointer=fopen(filename,option);
		if (filepointer != NULL)												//Checks if the file is found
		{
			if (fwrite(writeString,sizeof(char),strlen(writeString),filepointer) == strlen(writeString))		//Writes the string
				status = 1;														//Sets the status to 1
			fclose(filepointer);												//Closes the file
		}
	}
	return(status);																//Returns the status of the program
}


//############################################################################# READ A TEXT FILE
/**Routine for providing a text file stream.

\author JIC Vermaak*/
char* CHI_LIB::TextFileRead(char *fileName)
{
	FILE		*filePointer;													//Pointer to the file in memory
	char		*content = NULL;												//Pointer to the content being read
	//int			error;															//Error code
	int			count=0;														//Variable for getting file size

	if (fileName != NULL)														//Checks if a file name is indicated
	{
		//error = fopen_s(&filePointer,fileName,"rt");							//Opens the file for reading as a "TEXTFILE"
		filePointer=fopen(fileName,"rt");

		if (filePointer != NULL)												//If the file has successfully been opened
		{
			fseek(filePointer, 0, SEEK_END);									//Seeks the end of the file
			count = ftell(filePointer);											//Get the position of the position indicator
			rewind(filePointer);												//Puts the position indicator to the beginning of the file

			if (count > 0)														//If the file has some contents
			{
				content = (char *)malloc(sizeof(char) * (count+1));				//Allocates memory for the size of the content
				count = fread(content,sizeof(char),count,filePointer);			//Reads the file to content
				content[count] = '\0';											//Adds the null descriptor at the end of the file
			}
			fclose(filePointer);												//Closes the file
		}
		else
		{
			printf("Shader file not found\n");
		}

	}
	return content;																//Returns the content
}


