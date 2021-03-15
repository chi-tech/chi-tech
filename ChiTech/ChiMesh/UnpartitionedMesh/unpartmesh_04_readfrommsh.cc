#include "chi_unpartitioned_mesh.h"

#include "chi_log.h"
extern ChiLog& chi_log;
#include "chi_mpi.h"

//###################################################################
/**Reads an unpartitioned mesh from a wavefront .obj file.*/
void chi_mesh::UnpartitionedMesh::ReadFromMsh(const Options &options)
{
  //===================================================== Opening the file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Failed to open file: "<< options.file_name<<" in call "
      << "to ReadFromMsh \n";
    exit(EXIT_FAILURE);
  }

  chi_log.Log() << "Making Unpartitioned mesh from msh format file "
                << options.file_name;
  MPI_Barrier(MPI_COMM_WORLD);

  //===================================================== Reading every line and determining size
  std::string file_line;
  std::istringstream iss;
  const std::string node_section_name="$Nodes";
  const std::string elements_section_name = "$Elements";
  const std::string format_section_name = "$MeshFormat";

  //=================================================== Check the format of this input
  while (std::getline(file, file_line))
  {
    if ( format_section_name.compare(file_line)==0 )
      break;
  }
  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  double format;
  if ( !(iss >> format) )
  {
    chi_log.Log(LOG_ALLERROR)<<"Failed to read the file format.\n";
    exit(EXIT_FAILURE);
  }
  else if (format != 2.2)
  {
    chi_log.Log(LOG_ALLERROR)<<"Currently, only msh format 2.2 is supported.\n";
    exit(EXIT_FAILURE);
  }

  //=================================================== Find section with node information
  //                                                    and then read information
  file.seekg(0);
  while (std::getline(file, file_line))
  {
    if ( node_section_name.compare(file_line)==0 )
      break;
  }

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  int num_nodes;
  if ( !(iss >> num_nodes) )
  {
    chi_log.Log(LOG_ALLERROR)<<"Failed while trying to read the number of nodes.\n";
    exit(EXIT_FAILURE);
  }

  vertices.clear();
  vertices.resize(num_nodes);
  for (unsigned int ii = 0; ii < num_nodes; ++ii)
    vertices[ii] = new chi_mesh::Vertex;

  for (int n=0; n<num_nodes; n++)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    int vert_index;
    if ( !(iss >> vert_index) )
    {
      chi_log.Log(LOG_ALLERROR)<<"Failed to read vertex index.\n";
      exit(EXIT_FAILURE);
    }

    if (!(iss >> vertices[vert_index-1]->x >> vertices[vert_index-1]->y >> vertices[vert_index-1]->z))
    {
      chi_log.Log(LOG_ALLERROR)<<"Failed while reading the vertex coordinates.\n";
      exit(EXIT_FAILURE);
    }
  }

  //=================================================== Find the element listing section
  //                                                    and first read the boundary data
  file.seekg(0);
  while (std::getline(file, file_line))
  {
    if ( elements_section_name.compare(file_line)==0 )
      break;
  }

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  int num_elems;
  if (!(iss >> num_elems))
  {
    chi_log.Log(LOG_ALLERROR)<<"Failed to read number of elements.\n";
    exit(EXIT_FAILURE);
  }

  for (int n=0; n<num_elems; n++)
  {
    int elem_type, num_tags, physical_reg, tag, element_index;

    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    if ( !(iss >> element_index >> elem_type >> num_tags) )
    {
      chi_log.Log(LOG_ALLERROR)<<"Failed while reading element index, element type, and number of tags.\n";
      exit(EXIT_FAILURE);
    }

    if( !(iss>>physical_reg) )
    {
      chi_log.Log(LOG_ALLERROR)<<"Failed while reading physical region.\n";
      exit(EXIT_FAILURE);
    }

    for (int i=1; i<num_tags; i++)
      if( !(iss >> tag) )
      {
        chi_log.Log(LOG_ALLERROR)<<"Failed when reading tags.\n";
        exit(EXIT_FAILURE);
      }

    if (elem_type == 2)
    {
      raw_cells.push_back(new LightWeightCell(CellType::POLYGON));

      const int num_nodes = 3;

      int nodes[num_nodes];
      for (int i=0; i<num_nodes; i++)
        if ( !(iss >> nodes[i]) )
        {
          chi_log.Log(LOG_ALLERROR)<<"Failed when reading element node index.\n";
          exit(EXIT_FAILURE);
        }

      raw_cells.back()->vertex_ids.resize(num_nodes);
      for (int i=0; i<num_nodes; i++)
        raw_cells.back()->vertex_ids[i] = nodes[i]-1;

    }
    else if (elem_type == 3)
    {
      raw_cells.push_back(new LightWeightCell(CellType::POLYGON));

      const int num_nodes = 4;

      int nodes[num_nodes];
      for (int i=0; i<num_nodes; i++)
        if ( !(iss >> nodes[i]) )
        {
          chi_log.Log(LOG_ALLERROR)<<"Failed when reading element node index.\n";
          exit(EXIT_FAILURE);
        }

      raw_cells.back()->vertex_ids.resize(num_nodes);
      for (int i=0; i<num_nodes; i++)
        raw_cells.back()->vertex_ids[i] = nodes[i]-1;

    }
    else
    {
      continue;
    }

    const int total_nodes = raw_cells.back()->vertex_ids.size();

    for (int e=0; e<total_nodes; e++)
    {
      LightWeightFace face;

      face.vertex_ids.resize(2);

      face.vertex_ids[0] = raw_cells.back()->vertex_ids[e];

      if (e<total_nodes-1)
        face.vertex_ids[1] = raw_cells.back()->vertex_ids[e+1];
      else
        face.vertex_ids[1] = raw_cells.back()->vertex_ids[0];

      raw_cells.back()->faces.push_back(std::move(face));
    }
    raw_cells.back()->material_id = physical_reg;

  }

  file.close();

  //======================================== Always do this
  BuildMeshConnectivity();
  ComputeCentroids();
}



