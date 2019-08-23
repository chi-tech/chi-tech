#include "chi_lib.h"

#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#ifdef WINDOWS_ENV
    #include <afxres.h>
#endif


//#############################################################################
//####################### CREATE DIRECTORY ####################################
/**Creates a directory. Due to many anomolies with the WinAPI and directory
creation and deletion, this routine has been developed as a substitute routine
to enable multiple attempts as well as verification.

\param sPath				Folder path.
\param verifyCompletion		Flag indicating whether the routine should verify
							the existance of the folder. [Default: false].

\author JIC Vermaak*/
bool CHI_LIB::CreateDirectoryS(char* sPath,bool verifyCompletion)
{
#ifdef WINDOWS_ENV
	//===================================================== Attempt folder creation
	if (!CreateDirectory(sPath,NULL))
	{
		if (GetLastError() == ERROR_ALREADY_EXISTS) { /*DO NOTHING GO TO VERIFICATION*/}
		else
		{
			//============================================= Retry multiple times
			bool iterationFailed=true;
			for (int k=0; k<100; k++)
			{
				if (!CreateDirectory(sPath,NULL))
				{
					if (GetLastError() == ERROR_ALREADY_EXISTS) { iterationFailed=false; break; }
				}
				Sleep(5);
			}
			//============================================= If retries failed, post message
			if (iterationFailed)
			{
				char temp[1000];
				CHI_LIB::StringCopy(temp,"Error creating folder:\n");
				CHI_LIB::StringCatEx(temp,sPath);

				MessageBox(NULL,temp,"ERROR",MB_OK);
				return false;
			}
		}
	}

	//===================================================== Verify folder creation
	int k=0;
	if (verifyCompletion)
	{
		while ((GetFileAttributes(sPath) == INVALID_FILE_ATTRIBUTES) && (k<10000))
		{
			k++;
			Sleep(1);
		}
		if (k>=10000) {return false;}
	}
	return true;
#endif
  return false;
}







//#############################################################################
//####################### DELETE DIRECTORY ####################################
/**Recursively deletes the target folder and its contents. This routine is an
adaptation of the routine demonstrated at http://www.codeproject.com/KB/files/deletedir.aspx.

\param	sPath	The directory to be deleted.

\return Indicates whether or not the routine succeeded.

\author JIC Vermaak*/
bool CHI_LIB::DeleteDirectory(char* sPath)
{
#ifdef WINDOWS_ENV
    HANDLE hFind;  // file handle

    WIN32_FIND_DATA FindFileData;

    char DirPath[MAX_PATH];
    char FileName[MAX_PATH];

    CHI_LIB::StringCopy(DirPath,sPath);
    CHI_LIB::StringCatEx(DirPath,"\\*");    // searching all files

    CHI_LIB::StringCopy(FileName,sPath);
    CHI_LIB::StringCatEx(FileName,"\\");

    hFind = FindFirstFile(DirPath,&FindFileData); // find the first file
    if(hFind == INVALID_HANDLE_VALUE) return FALSE;
	CHI_LIB::StringCopy(DirPath,FileName);

    bool bSearch = true;
    while(bSearch)
	{ // until we finds an entry
        if(FindNextFile(hFind,&FindFileData))
		{
            //if(IsDots(FindFileData.cFileName)) continue;
			if ((CHI_LIB::StringCompare(FindFileData.cFileName,".") || CHI_LIB::StringCompare(FindFileData.cFileName,".."))) continue;

			CHI_LIB::StringCatEx(FileName,FindFileData.cFileName);

            if((FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
			{                // we have found a directory, recurse
                if(!DeleteDirectory(FileName))
				{
                    FindClose(hFind);
                    return FALSE; // directory couldn't be deleted
                }
                RemoveDirectory(FileName); // remove the empty directory
				CHI_LIB::StringCopy(FileName,DirPath);
            }
            else
			{
                if(FindFileData.dwFileAttributes & FILE_ATTRIBUTE_READONLY)
                    _chmod(FileName, _S_IWRITE); // change read-only file mode
                if(!DeleteFile(FileName))
				{  // delete the file
                    FindClose(hFind);
                    return FALSE;
                }
				CHI_LIB::StringCopy(FileName,DirPath);
            }
        }
        else
		{
            if(GetLastError() == ERROR_NO_MORE_FILES) // no more files there
            bSearch = false;
            else
			{
                // some error occured, close the handle and return FALSE
                FindClose(hFind);
                return FALSE;
            }

        }

    }
    FindClose(hFind);  // closing file handle

    return RemoveDirectory(sPath); // remove the empty directory
#endif
	return false;
}


//#############################################################################
//####################### DELETE DIRECTORY ####################################
/**Recursively deletes the contents of the folder. This routine is modified from
SmiLib::DeleteDirectory

\param	sPath	The directory to be cleared.

\return Indicates whether or not the routine succeeded.

\author JIC Vermaak*/
bool CHI_LIB::DeleteDirectoryContents(char* sPath)
{
#ifdef WINDOWS_ENV
    HANDLE hFind;  // file handle

    WIN32_FIND_DATA FindFileData;

    char DirPath[MAX_PATH];
    char FileName[MAX_PATH];

    CHI_LIB::StringCopy(DirPath,sPath);
    CHI_LIB::StringCatEx(DirPath,"\\*");    // searching all files

    CHI_LIB::StringCopy(FileName,sPath);
    CHI_LIB::StringCatEx(FileName,"\\");

    hFind = FindFirstFile(DirPath,&FindFileData); // find the first file
    if(hFind == INVALID_HANDLE_VALUE) return FALSE;
	CHI_LIB::StringCopy(DirPath,FileName);

    bool bSearch = true;
    while(bSearch)
	{ // until we finds an entry
        if(FindNextFile(hFind,&FindFileData))
		{
            //if(IsDots(FindFileData.cFileName)) continue;
			if ((CHI_LIB::StringCompare(FindFileData.cFileName,".") || CHI_LIB::StringCompare(FindFileData.cFileName,".."))) continue;

			CHI_LIB::StringCatEx(FileName,FindFileData.cFileName);

            if((FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
			{                // we have found a directory, recurse
                if(!DeleteDirectory(FileName))
				{
                    FindClose(hFind);
                    return FALSE; // directory couldn't be deleted
                }
                RemoveDirectory(FileName); // remove the empty directory
				CHI_LIB::StringCopy(FileName,DirPath);
            }
            else
			{
                if(FindFileData.dwFileAttributes & FILE_ATTRIBUTE_READONLY)
                    _chmod(FileName, _S_IWRITE); // change read-only file mode
                if(!DeleteFile(FileName))
				{  // delete the file
                    FindClose(hFind);
                    return FALSE;
                }
				CHI_LIB::StringCopy(FileName,DirPath);
            }
        }
        else
		{
            if(GetLastError() == ERROR_NO_MORE_FILES) // no more files there
            bSearch = false;
            else
			{
                // some error occured, close the handle and return FALSE
                FindClose(hFind);
                return FALSE;
            }

        }

    }
    FindClose(hFind);  // closing file handle

    return TRUE/*RemoveDirectory(sPath)*/; // remove the empty directory
#endif
	return false;
}




//#############################################################################
//###################### GET LIST OF FILES ####################################
/**Obtains a list of files, or the number of files in a folder. This routine
determines the amount of files in a folder if a fileList-pointer array is not
provided and returns the list of files if a fileList-pointer array is indeed
provided.

\param sPath	Name of the directory for which the filelist needs to be compiled.
\param fileList Pointer to an array of pointers for returning the file list.

\return			The number of files in this directory.

\author			JIC Vermaak*/
int CHI_LIB::GetListOfFiles(char* sPath,char* fileList[])
{
#ifdef WINDOWS_ENV
	int fileCount=0;

	HANDLE hFind;  // file handle

    WIN32_FIND_DATA FindFileData;

    char DirPath[MAX_PATH];
    char FileName[MAX_PATH];

    CHI_LIB::StringCopy(DirPath,sPath);
    CHI_LIB::StringCatEx(DirPath,"\\*");    // searching all files

    CHI_LIB::StringCopy(FileName,sPath);
    CHI_LIB::StringCatEx(FileName,"\\");

    hFind = FindFirstFile(DirPath,&FindFileData); // find the first file
    if(hFind == INVALID_HANDLE_VALUE) return 0;
	CHI_LIB::StringCopy(DirPath,FileName);

    bool bSearch = true;
    while(bSearch)
	{ // until we finds an entry
        if(FindNextFile(hFind,&FindFileData))
		{
            //if(IsDots(FindFileData.cFileName)) continue;
			if ((CHI_LIB::StringCompare(FindFileData.cFileName,".") || CHI_LIB::StringCompare(FindFileData.cFileName,".."))) continue;

			CHI_LIB::StringCopy(FileName,FindFileData.cFileName);

            if((FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
			{
				// we have found a directory, skip
            }
            else
			{
				fileCount++;
				if (fileList != NULL)
				{
					fileList[fileCount-1] = new char[1000];
					CHI_LIB::StringCopy(fileList[fileCount-1],FindFileData.cFileName);
				}
            }
        }
        else
		{
            if(GetLastError() == ERROR_NO_MORE_FILES) // no more files there
            bSearch = false;
            else
			{
                // some error occured, close the handle and return FALSE
                FindClose(hFind);
                return fileCount;
            }

        }

    }
    FindClose(hFind);  // closing file handle

    return fileCount; // remove the empty directory
#endif
	return false;
}
