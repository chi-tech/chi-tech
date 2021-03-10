#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "chi_log.h"
#include "chi_mpi.h"
extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

#include <sys/stat.h>
#include <fstream>
#include <cstring>

//###################################################################
/**Writes phi_old to restart file.*/
void LinearBoltzmann::Solver::WriteRestartData(std::string folder_name,
                                               std::string file_base)
{
  typedef struct stat Stat;
  Stat st;

  //======================================== Make sure folder exists
  if (chi_mpi.location_id == 0)
  {
    if (stat(folder_name.c_str(),&st) != 0) //if not exist, make it
      if ( (mkdir(folder_name.c_str(),S_IRWXU | S_IRWXG | S_IRWXO) != 0) and
           (errno != EEXIST) )
      {
        chi_log.Log(LOG_0WARNING)
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
  sprintf(location_cstr,"%d.r",chi_mpi.location_id);

  std::string file_name = folder_name + std::string("/") +
                          file_base + std::string(location_cstr);

  std::ofstream ofile;
  ofile.open(file_name, std::ios::out | std::ios::binary | std::ios::trunc);

  if (not ofile.is_open())
  {
    chi_log.Log(LOG_ALLERROR)
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
    chi_log.Log(LOG_0)
      << "Successfully wrote restart data: "
      << folder_name + std::string("/") +
         file_base + std::string("X.r");
  else
    chi_log.Log(LOG_0ERROR)
      << "Failed to write restart data: "
      << folder_name + std::string("/") +
         file_base + std::string("X.r");
}

//###################################################################
/**Read phi_old from restart file.*/
void LinearBoltzmann::Solver::ReadRestartData(std::string folder_name,
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
  sprintf(location_cstr,"%d.r",chi_mpi.location_id);

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
    chi_log.Log(LOG_0) << "Successfully read restart data";
  else
    chi_log.Log(LOG_0ERROR)
      << "Failed to read restart data: "
      << folder_name + std::string("/") +
         file_base + std::string("X.r");
}

//###################################################################
/**Prints the groupset's angular fluxes to file.*/
void LinearBoltzmann::Solver::
  WriteGroupsetAngularFluxes(const LBSGroupset& groupset,
                             const std::string& file_base)
{
  std::string file_name =
    file_base + std::to_string(chi_mpi.location_id) + ".data";

  //============================================= Open file
  std::ofstream file(file_name,
                     std::ofstream::binary | //binary file
                     std::ofstream::out |    //no accidental reading
                     std::ofstream::trunc);  //clear file contents when opened

  //============================================= Check file is open
  if (not file.is_open())
  {
    chi_log.Log(LOG_ALLWARNING)
      << __FUNCTION__ << "Failed to open " << file_name;
    return;
  }

  //============================================= Write header
  std::string header_info =
    "Chi-Tech LinearBoltzmann::Groupset angular flux file\n"
    "Header size: 320 bytes\n"
    "Structure(type-info):\n"
    "size_t-num_local_nodes\n"
    "size_t-num_angles\n"
    "size_t-num_groups\n"
    "size_t-num_records\n"
    "Each record:\n"
    "size_t-cell_global_id\n"
    "unsigned int-node_number\n"
    "unsigned int-angle_num\n"
    "unsigned int-group_num\n"
    "double-angular_flux\n";

  int header_size = (int)header_info.length();

  char header_bytes[320];
  memset(header_bytes, '-', 320);
  strncpy(header_bytes, header_info.c_str(),std::min(header_size,319));
  header_bytes[319]='\0';

  file << header_bytes;

  auto NODES_ONLY = ChiMath::UNITARY_UNKNOWN_MANAGER;

  //============================================= Get relevant items
  auto fe = std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(discretization);
  if (not fe) {file.close(); return;}
  size_t num_local_nodes = discretization->GetNumLocalDOFs(grid,NODES_ONLY);
  size_t num_angles      = groupset.quadrature->abscissae.size();
  size_t num_groups      = groupset.groups.size();
  size_t num_local_dofs  = groupset.num_psi_unknowns_local;
  auto   dof_handler     = groupset.psi_uk_man;

  //============================================= Write quantities
  file.write((char*)&num_local_nodes,sizeof(size_t));
  file.write((char*)&num_angles     ,sizeof(size_t));
  file.write((char*)&num_groups     ,sizeof(size_t));
  file.write((char*)&num_local_dofs ,sizeof(size_t));

  //============================================= Write per node data
  size_t dof_count=0;
  for (const auto& cell : grid->local_cells)
  {
    const auto cell_fe_mapping = fe->GetCellMappingFE(cell.local_id);
    for (unsigned int i=0; i<cell_fe_mapping->num_nodes; ++i)
      for (unsigned int n=0; n<num_angles; ++n)
        for (unsigned int g=0; g<num_groups; ++g)
        {
          if (++dof_count > num_local_dofs) goto close_file;
          uint64_t dof_map = fe->MapDOFLocal(cell,i,dof_handler,n,g);
          double value = groupset.psi_new_local[dof_map];

          file.write((char*)&cell.global_id,sizeof(size_t));
          file.write((char*)&i             ,sizeof(unsigned int));
          file.write((char*)&n             ,sizeof(unsigned int));
          file.write((char*)&g             ,sizeof(unsigned int));
          file.write((char*)&value         ,sizeof(double));
        }
  }

  //============================================= Clean-up
close_file:
  file.close();
}

//###################################################################
/**Prints the groupset's angular fluxes to file.*/
void LinearBoltzmann::Solver::
  ReadGroupsetAngularFluxes(const LBSGroupset& groupset,
                           const std::string& file_base)
{
  std::string file_name =
    file_base + std::to_string(chi_mpi.location_id) + ".data";

  //============================================= Open file
  std::ifstream file(file_name,
                     std::ofstream::binary | //binary file
                     std::ofstream::in);     //no accidental writing

  //============================================= Check file is open
  if (not file.is_open())
  {
    chi_log.Log(LOG_ALLWARNING)
      << __FUNCTION__ << "Failed to open " << file_name;
    return;
  }

  size_t num_local_nodes;
  size_t num_angles     ;
  size_t num_groups     ;
  size_t num_local_dofs ;

  std::vector<double> psi;

  //============================================= Read header
  char header_bytes[320]; header_bytes[319] = '\0';
  file.read(header_bytes,319);

  file.read((char*)&num_local_nodes,sizeof(size_t));
  file.read((char*)&num_angles     ,sizeof(size_t));
  file.read((char*)&num_groups     ,sizeof(size_t));
  file.read((char*)&num_local_dofs ,sizeof(size_t));

  chi_log.Log(LOG_ALL) << header_bytes;

  {
    std::stringstream outstr;
    outstr << "num_local_nodes: " << num_local_nodes << "\n";
    outstr << "num_angles     : " << num_angles      << "\n";
    outstr << "num_groups     : " << num_groups      << "\n";
    outstr << "num_local_dofs : " << num_local_dofs  << "\n";
    chi_log.Log(LOG_ALL) << outstr.str();
  }

  psi.reserve(num_local_dofs);
  std::set<uint64_t> cells_touched;
  for (size_t dof=0; dof<num_local_dofs; ++dof)
  {
    uint64_t     cell_global_id;
    unsigned int node;
    unsigned int angle_num;
    unsigned int group;
    double       psi_value;

    file.read((char*)&cell_global_id,sizeof(uint64_t));
    file.read((char*)&node          ,sizeof(unsigned int));
    file.read((char*)&angle_num     ,sizeof(unsigned int));
    file.read((char*)&group         ,sizeof(unsigned int));
    file.read((char*)&psi_value     ,sizeof(double));

    cells_touched.insert(cell_global_id);
  }

  chi_log.Log(LOG_ALL) << "Number of cells read: " << cells_touched.size();

  //============================================= Clean-up
  file.close();
}