#include "chi_log.h"
#include <chi_mpi.h>


extern CHI_MPI chi_mpi;

//###################################################################
/** Default constructor*/
CHI_LOG::CHI_LOG()
{
  verbosity = LOG_0VERBOSE_0;

}

//###################################################################
/** Makes a log entry.*/
CHI_LOG_STREAM CHI_LOG::Log(LOG_LVL level)
{


  switch (level)
  {
    case LOG_0:
    {
      if (chi_mpi.location_id == 0)
      {
        std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
        return CHI_LOG_STREAM(&std::cout, header);
      }
      else
      {
        std::string header = " ";
        return CHI_LOG_STREAM(&dummy_stream, header);
      }
    }
    case LOG_0WARNING:
    {
      if (chi_mpi.location_id == 0)
      {
        std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
        header += "**WARNING** ";
        return CHI_LOG_STREAM(&std::cout, header);
      }
      else
      {
        std::string header = " ";
        return CHI_LOG_STREAM(&dummy_stream, header);
      }
    }
    case LOG_0ERROR:
    {
      if (chi_mpi.location_id == 0)
      {
        std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
        header += "**!**ERROR**!** ";
        return CHI_LOG_STREAM(&std::cerr, header);
      }
      else
      {
        std::string header = " ";
        return CHI_LOG_STREAM(&dummy_stream, header);
      }
    }
    case LOG_0VERBOSE_0:
    case LOG_0VERBOSE_1:
    case LOG_0VERBOSE_2:
    {
      if ((chi_mpi.location_id == 0) && (verbosity >= level))
      {
        std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
        return CHI_LOG_STREAM(&std::cout, header);
      }
      else
      {
        std::string header = " ";
        return CHI_LOG_STREAM(&dummy_stream, header);
      }
    }
    case LOG_ALL:
    {
      std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
      return CHI_LOG_STREAM(&std::cout, header);
    }
    case LOG_ALLWARNING:
    {
      std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
      header += "**WARNING** ";
      return CHI_LOG_STREAM(&std::cout, header);
    }
    case LOG_ALLERROR:
    {
      std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
      header += "**!**ERROR**!** ";
      return CHI_LOG_STREAM(&std::cerr, header);
    }

    case LOG_ALLVERBOSE_0:
    case LOG_ALLVERBOSE_1:
    case LOG_ALLVERBOSE_2:
    {
      if (verbosity >= (level-6))
      {
        std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
        return CHI_LOG_STREAM(&std::cout, header);
      }
      else
      {
        std::string header = " ";
        return CHI_LOG_STREAM(&dummy_stream, header);
      }
    }
  }


  std::string header = " ";
  return CHI_LOG_STREAM(&dummy_stream, header);
}


//###################################################################
/** Sets the verbosity level.*/
void CHI_LOG::SetVerbosity(int int_level)
{
  if (int_level == 0)
  {
    verbosity = LOG_0VERBOSE_0;
  }
  else if (int_level == 1)
  {
    verbosity = LOG_0VERBOSE_1;
  }
  else if (int_level == 2)
  {
    verbosity = LOG_0VERBOSE_2;
  }
}

//###################################################################
/** Gets the current verbosity level.*/
int CHI_LOG::GetVerbosity()
{
  return verbosity;
}
