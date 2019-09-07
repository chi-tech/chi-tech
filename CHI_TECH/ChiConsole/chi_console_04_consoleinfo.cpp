#include "chi_console.h"


#include <iostream>
#include <cstdio>
#include <sstream>

#ifdef WINDOWS_ENV
#include <windows.h>
#endif

//############################################################################# Read console size
/** Gets the amount of characters in the current buffer.
\author Jan*/
int CHI_CONSOLE::GetNumCharsInConsoleBuffer()
{
#ifdef WINDOWS_ENV
	CONSOLE_SCREEN_BUFFER_INFO buffer_info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &buffer_info);
	this->numberOfLines=buffer_info.dwCursorPosition.Y;
	this->xSize=buffer_info.dwSize.X;
	return (int)((buffer_info.dwSize.X * ( buffer_info.dwCursorPosition.Y + 1)) - (buffer_info.dwSize.X - (buffer_info.dwCursorPosition.X  + 1)));

#endif
	return 0;
}





//############################################################################# Provides a copy of the console
/**
\author Jan*/
void CHI_CONSOLE::CopyConsole(char* destination, int lineNumber,int xSize)
{
#ifdef WINDOWS_ENV
	DWORD num_character_read = 0;
	COORD first_char_to_read = {0,this->numberOfLines+lineNumber};
	ReadConsoleOutputCharacter(GetStdHandle(STD_OUTPUT_HANDLE), destination, xSize, first_char_to_read, &num_character_read);
	//this->buffer[2000-1] = '\0';
#endif
}
