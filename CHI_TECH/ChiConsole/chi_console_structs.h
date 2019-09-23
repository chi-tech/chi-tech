#ifndef CHI_CONSOLE_STRUCTS_H
#define CHI_CONSOLE_STRUCTS_H


//=============================================================================
/**Stores data relevant to events*/
struct CSTEvent
{
  char  eventTitle[100];
  char  sPar[6][100];
  int   iPar[6];
  float fPar[6];
  bool  bPar[6];

  CSTEvent()
  {
    eventTitle[0] = '\0';
    for (int k = 0; k < 6; ++k)
    {
      sPar[k][0] = '\0';
      iPar[k]    = 0;
      fPar[k]    = 0.0;
      bPar[k]    = false;
    }
  }
};

//=============================================================================
/**Simple structure for memory usage.*/
struct CSTMemory
{
  double memory_bytes;
  double memory_kbytes;
  double memory_mbytes;
  double memory_gbytes;

  CSTMemory()
  {
    memory_bytes  = 0.0;
    memory_kbytes = 0.0;
    memory_mbytes = 0.0;
    memory_gbytes = 0.0;
  }

  explicit CSTMemory(double in_mem)
  {
    memory_bytes  = in_mem;
    memory_kbytes = in_mem/1024.0;
    memory_mbytes = in_mem/1024.0/1024.0;
    memory_gbytes = in_mem/1024.0/1024.0;
  }

  CSTMemory& operator=(const CSTMemory& in_struct)
  {
    memory_bytes  = in_struct.memory_bytes ;
    memory_kbytes = in_struct.memory_kbytes;
    memory_mbytes = in_struct.memory_mbytes;
    memory_gbytes = in_struct.memory_gbytes;
    return *this;
  }
};

#endif
