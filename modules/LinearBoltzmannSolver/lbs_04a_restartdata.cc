#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"
;


#include <sys/stat.h>
#include <fstream>
#include <cstring>

//###################################################################
/**Writes phi_old to restart file.*/
void lbs::SteadySolver::WriteRestartData(std::string folder_name,
                                         std::string file_base)
{
  typedef struct stat Stat;
  Stat st;

  //======================================== Make sure folder exists
  if (chi::mpi.location_id == 0)
  {
    if (stat(folder_name.c_str(),&st) != 0) //if not exist, make it
      if ( (mkdir(folder_name.c_str(),S_IRWXU | S_IRWXG | S_IRWXO) != 0) and
           (errno != EEXIST) )
      {
        chi::log.Log0Warning()
          << "Failed to create restart directory: " << folder_name;
        return;
      }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //======================================== Create files
  //This step might fail for specific locations and
  //can create quite a messy output if we print it all.
  //We also need to consolidate the error to determine if
  //the process as whole succeeded.
  bool location_succeeded = true;
  char location_cstr[20];
  snprintf(location_cstr,20,"%d.r",chi::mpi.location_id);

  std::string file_name = folder_name + std::string("/") +
                          file_base + std::string(location_cstr);

  std::ofstream ofile;
  ofile.open(file_name, std::ios::out | std::ios::binary | std::ios::trunc);

  if (not ofile.is_open())
  {
    chi::log.LogAllError()
      << "Failed to create restart file: " << file_name;
    ofile.close();
    location_succeeded = false;
  }
  else
  {
    size_t phi_old_size = phi_old_local.size();
    ofile.write((char*)&phi_old_size, sizeof(size_t));
    for (auto val : phi_old_local)
      ofile.write((char*)&val, sizeof(double));

    ofile.close();
  }

  //======================================== Wait for all processes
  //                                         then check success status
  MPI_Barrier(MPI_COMM_WORLD);
  bool global_succeeded = true;
  MPI_Allreduce(&location_succeeded,   //Send buffer
                &global_succeeded,     //Recv buffer
                1,                     //count
                MPI_CXX_BOOL,          //Data type
                MPI_LAND,              //Operation - Logical and
                MPI_COMM_WORLD);       //Communicator

  //======================================== Write status message
  if (global_succeeded)
    chi::log.Log()
      << "Successfully wrote restart data: "
      << folder_name + std::string("/") +
         file_base + std::string("X.r");
  else
    chi::log.Log0Error()
      << "Failed to write restart data: "
      << folder_name + std::string("/") +
         file_base + std::string("X.r");
}

//###################################################################
/**Read phi_old from restart file.*/
void lbs::SteadySolver::ReadRestartData(std::string folder_name,
                                        std::string file_base)
{
  MPI_Barrier(MPI_COMM_WORLD);

  //======================================== Open files
  //This step might fail for specific locations and
  //can create quite a messy output if we print it all.
  //We also need to consolidate the error to determine if
  //the process as whole succeeded.
  bool location_succeeded = true;
  char location_cstr[20];
  snprintf(location_cstr,20,"%d.r",chi::mpi.location_id);

  std::string file_name = folder_name + std::string("/") +
                          file_base + std::string(location_cstr);

  std::ifstream ifile;
  ifile.open(file_name, std::ios::in | std::ios::binary );

  if (not ifile.is_open())
  {
    ifile.close();
    location_succeeded = false;
  }
  else
  {
    size_t number_of_unknowns;
    ifile.read((char*)&number_of_unknowns, sizeof(size_t));

    if (number_of_unknowns != phi_old_local.size())
    {
      location_succeeded = false;
      ifile.close();
    }
    else
    {
      std::vector<double> temp_phi_old(phi_old_local.size(),0.0);

      size_t v=0;
      while (not ifile.eof())
      {
        ifile.read((char*)&temp_phi_old[v], sizeof(double));
        ++v;
      }

      if (v != (number_of_unknowns+1))
      {
        location_succeeded = false;
        ifile.close();
      }
      else
        phi_old_local = std::move(temp_phi_old);

      ifile.close();
    }
  }

  //======================================== Wait for all processes
  //                                         then check success status
  MPI_Barrier(MPI_COMM_WORLD);
  bool global_succeeded = true;
  MPI_Allreduce(&location_succeeded,   //Send buffer
                &global_succeeded,     //Recv buffer
                1,                     //count
                MPI_CXX_BOOL,          //Data type
                MPI_LAND,              //Operation - Logical and
                MPI_COMM_WORLD);       //Communicator

  //======================================== Write status message
  if (global_succeeded)
    chi::log.Log() << "Successfully read restart data";
  else
    chi::log.Log0Error()
      << "Failed to read restart data: "
      << folder_name + std::string("/") +
         file_base + std::string("X.r");
}
