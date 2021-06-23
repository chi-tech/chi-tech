#include "lbs_linear_boltzmann_solver.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include <fstream>
#include <cstring>

//###################################################################
/**Makes a source-moments vector from scattering and fission based
 * on the latest phi-solution.*/
void LinearBoltzmann::Solver::
  MakeSourceMomentsFromPhi(std::vector<double> &source_moments)
{
  size_t num_local_dofs = discretization->GetNumLocalDOFs(flux_moments_uk_man);

  source_moments.assign(num_local_dofs,0.0);
  for (auto& groupset : group_sets)
  {
    SetSource(groupset, APPLY_AGS_SCATTER_SOURCE | APPLY_WGS_SCATTER_SOURCE |
                        APPLY_FISSION_SOURCE);
    ScopedCopySTLvectors(groupset, q_moments_local, source_moments);
  }
}


//###################################################################
/**Writes a given flux-moments vector to file.*/
void LinearBoltzmann::Solver::
  WriteFluxMoments(const std::string &file_base,
                   const std::vector<double>& flux_moments)
{
  std::string file_name =
    file_base + std::to_string(chi_mpi.location_id) + ".data";

  //============================================= Open file
  chi_log.Log() << "Writing flux-moments file " << file_name;
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
    "Chi-Tech LinearBoltzmann: Flux moments file\n"
    "Header size: 320 bytes\n"
    "Structure(type-info):\n"
    "size_t-num_local_nodes\n"
    "size_t-num_moments\n"
    "size_t-num_groups\n"
    "size_t-num_records\n"
    "Each record:\n"
    "size_t-cell_global_id\n"
    "unsigned int-node_number\n"
    "unsigned int-moment_num\n"
    "unsigned int-group_num\n"
    "double-flux_moment_value\n";

  int header_size = (int)header_info.length();

  char header_bytes[320];
  memset(header_bytes, '-', 320);
  strncpy(header_bytes, header_info.c_str(),std::min(header_size,319));
  header_bytes[319]='\0';

  file << header_bytes;

  //============================================= Get relevant items
  auto NODES_ONLY = ChiMath::UNITARY_UNKNOWN_MANAGER;
  auto& sdm = discretization;
  size_t num_local_nodes = discretization->GetNumLocalDOFs(NODES_ONLY);
  size_t num_moments_t   = static_cast<size_t>(num_moments);
  size_t num_groups      = groups.size();
  size_t num_local_dofs  = discretization->GetNumLocalDOFs(flux_moments_uk_man);

  //============================================= Write num_ quantities
  file.write((char*)&num_local_nodes,sizeof(size_t));
  file.write((char*)&num_moments_t  ,sizeof(size_t));
  file.write((char*)&num_groups     ,sizeof(size_t));
  file.write((char*)&num_local_dofs ,sizeof(size_t));

  //============================================= Write per dof data
  try {
    for (const auto& cell : grid->local_cells)
      for (unsigned int i=0; i<sdm->GetCellNumNodes(cell); ++i)
        for (unsigned int m=0; m<num_moments; ++m)
          for (unsigned int g=0; g<num_groups; ++g)
          {
            uint64_t dof_map = sdm->MapDOFLocal(cell,i,flux_moments_uk_man,m,g);
            double value = flux_moments.at(dof_map);

            file.write((char*)&cell.global_id,sizeof(size_t));
            file.write((char*)&i             ,sizeof(unsigned int));
            file.write((char*)&m             ,sizeof(unsigned int));
            file.write((char*)&g             ,sizeof(unsigned int));
            file.write((char*)&value         ,sizeof(double));
          }
  }
  catch (const std::out_of_range& e)
  {
    chi_log.Log(LOG_ALLWARNING) << __FUNCTION__ << ": The given flux_moments-"
                               << "vector was accessed out of range.";
  }

  //============================================= Clean-up
  file.close();
}


//###################################################################
/**Writes a flux-moments vector from a file in the specified vector.*/
void LinearBoltzmann::Solver::ReadFluxMoments(const std::string &file_base,
                                              std::vector<double>& flux_moments,
                                              bool single_file/*=false*/)
{
  std::string file_name =
    file_base + std::to_string(chi_mpi.location_id) + ".data";
  if (single_file)
    file_name = file_base + ".data";

  //============================================= Open file
  chi_log.Log() << "Reading flux-moments file " << file_name;
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

  //============================================= Get relevant items
  auto NODES_ONLY = ChiMath::UNITARY_UNKNOWN_MANAGER;
  auto& sdm = discretization;
  size_t num_local_nodes = discretization->GetNumLocalDOFs(NODES_ONLY);
  size_t num_moments_t   = static_cast<size_t>(num_moments);
  size_t num_groups      = groups.size();
  size_t num_local_dofs  = discretization->GetNumLocalDOFs(flux_moments_uk_man);

  size_t file_num_local_nodes;
  size_t file_num_moments    ;
  size_t file_num_groups     ;
  size_t file_num_local_dofs ;

  flux_moments.assign(num_local_dofs,0.0);

  //============================================= Read header
  char header_bytes[320]; header_bytes[319] = '\0';
  file.read(header_bytes,319);

  file.read((char*)&file_num_local_nodes, sizeof(size_t));
  file.read((char*)&file_num_moments    , sizeof(size_t));
  file.read((char*)&file_num_groups     , sizeof(size_t));
  file.read((char*)&file_num_local_dofs , sizeof(size_t));

  //============================================= Check compatibility
  if (not single_file)
    if (file_num_local_nodes != num_local_nodes or
        file_num_moments     != num_moments_t   or
        file_num_groups      != num_groups      or
        file_num_local_dofs  != num_local_dofs)
    {
      std::stringstream outstr;
      outstr << "num_local_nodes: " << file_num_local_nodes << " vs "
                                    << num_local_nodes << "\n";
      outstr << "num_moments    : " << file_num_moments << " vs "
                                    << num_moments_t << "\n";
      outstr << "num_groups     : " << file_num_groups << " vs "
                                    << num_groups << "\n";
      outstr << "num_local_dofs : " << file_num_local_dofs << " vs "
                                    << num_local_dofs << "\n";
      chi_log.Log(LOG_ALL)
        << "Incompatible DOF data found in file " << file_name << "\n"
        << outstr.str();
      file.close();
      return;
    }

  //============================================= Commit to reading the file
  std::set<uint64_t> cells_touched;
  for (size_t dof=0; dof < file_num_local_dofs; ++dof)
  {
    uint64_t     cell_global_id;
    unsigned int node;
    unsigned int moment;
    unsigned int group;
    double       flux_value;

    file.read((char*)&cell_global_id,sizeof(uint64_t));
    file.read((char*)&node          ,sizeof(unsigned int));
    file.read((char*)&moment        ,sizeof(unsigned int));
    file.read((char*)&group         ,sizeof(unsigned int));
    file.read((char*)&flux_value    ,sizeof(double));

    if (grid->IsCellLocal(cell_global_id))
    {
      cells_touched.insert(cell_global_id);

      const auto& cell = grid->cells[cell_global_id];

      size_t dof_map = sdm->MapDOFLocal(cell,node,flux_moments_uk_man,moment,group);

      flux_moments.at(dof_map) = flux_value;
    }
  }

  //============================================= Clean-up
  file.close();
}