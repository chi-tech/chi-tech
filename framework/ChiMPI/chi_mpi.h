#ifndef CHI_MPI_H
#define CHI_MPI_H

#include <mpi.h>
#include "../ChiMesh/chi_mesh.h"
#include "chi_runtime.h"

namespace chi_objects
{

  //################################################################### Class def
  /**An object for storing various MPI states.*/
  class MPI_Info
  {
  private:
    int location_id_ = 0;
    int process_count_ = 1;

    bool location_id_set_ = false;
    bool process_count_set_ = false;

  public:
    const int& location_id = location_id_;     ///< Current process rank.
    const int& process_count = process_count_; ///< Total number of processes.

  private:
    static MPI_Info instance;
    MPI_Info() = default;

  public:
    MPI_Info(const MPI_Info&) = delete;           //Deleted copy constructor
    MPI_Info operator=(const MPI_Info&) = delete; //Deleted assigment operator

  private:
    friend int chi::Initialize(int argc, char** argv);
    void SetLocationID(int in_location_id)
    {
      if (not location_id_set_)
        location_id_ = in_location_id;
      location_id_set_ = true;
    }
    void SetProcessCount(int in_process_count)
    {
      if (not process_count_set_)
        process_count_ = in_process_count;
      process_count_set_ = true;
    }

  public:
    static MPI_Info& GetInstance() noexcept {return instance;}
  };
}//namespace chi_objects



#endif