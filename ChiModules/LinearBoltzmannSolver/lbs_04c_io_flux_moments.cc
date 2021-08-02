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
  for (auto& groupset : groupsets)
  {
    SetSource(groupset,
              source_moments,
              APPLY_AGS_SCATTER_SOURCE | APPLY_WGS_SCATTER_SOURCE |
              APPLY_AGS_FISSION_SOURCE | APPLY_WGS_FISSION_SOURCE);
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
  chi_log.Log() << "Writing flux-moments to files with base-name " << file_base
                << " and extension .data";
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
    "Header size: 500 bytes\n"
    "Structure(type-info):\n"
    "size_t num_local_nodes\n"
    "size_t num_moments\n"
    "size_t num_groups\n"
    "size_t num_records\n"
    "size_t num_cells\n"
    "Each cell:\n"
    "  uint64_t cell_global_id\n"
    "  size_t   num_nodes\n"
    "  Each node:\n"
    "    double   x_position\n"
    "    double   y_position\n"
    "    double   z_position\n"
    "Each record:\n"
    "  uint64_t     cell_global_id\n"
    "  unsigned int node_number\n"
    "  unsigned int moment_num\n"
    "  unsigned int group_num\n"
    "  double       flux_moment_value\n";

  int header_size = (int)header_info.length();

  char header_bytes[500];
  memset(header_bytes, '-', 500);
  strncpy(header_bytes, header_info.c_str(),std::min(header_size,499));
  header_bytes[499]='\0';

  file << header_bytes;

  //============================================= Get relevant items
  auto NODES_ONLY = ChiMath::UNITARY_UNKNOWN_MANAGER;
  auto& sdm = discretization;
  size_t num_local_nodes = discretization->GetNumLocalDOFs(NODES_ONLY);
  size_t num_moments_t   = static_cast<size_t>(num_moments);
  size_t num_groups      = groups.size();
  size_t num_local_dofs  = discretization->GetNumLocalDOFs(flux_moments_uk_man);
  size_t num_local_cells = grid->local_cells.size();

  //============================================= Write num_ quantities
  file.write((char*)&num_local_nodes,sizeof(size_t));
  file.write((char*)&num_moments_t  ,sizeof(size_t));
  file.write((char*)&num_groups     ,sizeof(size_t));
  file.write((char*)&num_local_dofs ,sizeof(size_t));
  file.write((char*)&num_local_cells,sizeof(size_t));

  //============================================= Write nodal positions for
  //                                              each cell
  for (const auto& cell : grid->local_cells)
  {
    file.write((char *) &cell.global_id, sizeof(uint64_t));

    size_t num_nodes      = discretization->GetCellNumNodes(cell);
    file.write((char *) &num_nodes, sizeof(size_t));

    auto   node_locations = discretization->GetCellNodeLocations(cell);
    for (const auto& node : node_locations)
    {
      file.write((char *) &node.x, sizeof(double));
      file.write((char *) &node.y, sizeof(double));
      file.write((char *) &node.z, sizeof(double));
    }//for node
  }//for cell

  //============================================= Write per dof data
  try {
    for (const auto& cell : grid->local_cells)
      for (unsigned int i=0; i<sdm->GetCellNumNodes(cell); ++i)
        for (unsigned int m=0; m<num_moments_t; ++m)
          for (unsigned int g=0; g<num_groups; ++g)
          {
            uint64_t dof_map = sdm->MapDOFLocal(cell,i,flux_moments_uk_man,m,g);
            double value = flux_moments.at(dof_map);

            file.write((char*)&cell.global_id,sizeof(uint64_t));
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
  size_t num_local_cells = grid->local_cells.size();

  size_t file_num_local_nodes;
  size_t file_num_moments    ;
  size_t file_num_groups     ;
  size_t file_num_local_dofs ;
  size_t file_num_local_cells;

  flux_moments.assign(num_local_dofs,0.0);

  //============================================= Read header
  char header_bytes[500]; header_bytes[499] = '\0';
  file.read(header_bytes,499);

  file.read((char*)&file_num_local_nodes, sizeof(size_t));
  file.read((char*)&file_num_moments    , sizeof(size_t));
  file.read((char*)&file_num_groups     , sizeof(size_t));
  file.read((char*)&file_num_local_dofs , sizeof(size_t));
  file.read((char*)&file_num_local_cells, sizeof(size_t));

  //============================================= Check compatibility
  if (not single_file)
    if (file_num_local_nodes != num_local_nodes or
        file_num_moments     != num_moments_t   or
        file_num_groups      != num_groups      or
        file_num_local_dofs  != num_local_dofs  or
        file_num_local_cells != num_local_cells)
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
      outstr << "num_local_cells: " << file_num_local_cells << " vs "
                                    << num_local_cells << "\n";
      chi_log.Log(LOG_ALL)
        << "Incompatible DOF data found in file " << file_name << "\n"
        << "File data vs system:\n" << outstr.str();
      file.close();
      return;
    }

  //============================================= Read cell nodal locations
  std::map<uint64_t, std::map<size_t,size_t>> file_cell_nodal_mapping;
  for (size_t c=0; c < file_num_local_cells; ++c)
  {
    //============================ Read cell-id and num_nodes
    uint64_t cell_global_id;
    size_t   num_nodes;

    file.read((char*)&cell_global_id, sizeof(uint64_t));
    file.read((char*)&num_nodes     , sizeof(size_t));

    //============================ Read node locations
    std::vector<chi_mesh::Vector3> file_node_locations;
    file_node_locations.reserve(num_nodes);
    for (size_t n=0; n < num_nodes; ++n)
    {
      double x,y,z;
      file.read((char*)&x, sizeof(double));
      file.read((char*)&y, sizeof(double));
      file.read((char*)&z, sizeof(double));

      file_node_locations.emplace_back(x,y,z);
    }//for file node n

    if (not grid->IsCellLocal(cell_global_id)) continue;

    const auto& cell = grid->cells[cell_global_id];

    //================ Now map file nodes to system nodes
    auto system_node_locations = discretization->GetCellNodeLocations(cell);
    std::map<size_t,size_t> mapping;

    //Check num_nodes equal
    if (system_node_locations.size() != num_nodes)
      throw std::logic_error(std::string(__FUNCTION__) +
           ": Incompatible number of nodes for a cell was encountered. Mapping "
           "could not be performed.");

    bool mapping_successful = true; //Assume true, now try to disprove

    const auto &sys_nodes = system_node_locations;
    const auto &file_nodes = file_node_locations;

    for (size_t n = 0; n < num_nodes; ++n)
    {
      bool mapping_found = false;
      for (size_t m = 0; m < num_nodes; ++m)
        if ((sys_nodes[m] - file_nodes[n]).NormSquare() < 1.0e-12) {
          mapping[n] = m; mapping_found = true; }

      if (not mapping_found) {
        mapping_successful = false; break; }
    }//for n

    if (not mapping_successful)
      throw std::logic_error(std::string(__FUNCTION__) +
           ": Incompatible node locations for a cell was encountered. Mapping "
           "unsuccessful.");

    file_cell_nodal_mapping[cell_global_id] = std::move(mapping);
  }//for c (cell in file)

  //============================================= Commit to reading the file
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
      const auto& cell         = grid->cells[cell_global_id];
      const auto& node_mapping = file_cell_nodal_mapping.at(cell_global_id);

      size_t node_mapped = node_mapping.at(node);
//      auto node_mapped = node;

      size_t dof_map = sdm->MapDOFLocal(cell,
                                        node_mapped,
                                        flux_moments_uk_man,
                                        moment,
                                        group);

      flux_moments.at(dof_map) = flux_value;
    }//if cell is local
  }//for dof

  //============================================= Clean-up
  file.close();
}