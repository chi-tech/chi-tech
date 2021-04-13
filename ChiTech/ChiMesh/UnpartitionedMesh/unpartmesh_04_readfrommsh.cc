#include "chi_unpartitioned_mesh.h"

#include "chi_log.h"
extern ChiLog& chi_log;
#include "chi_mpi.h"

//###################################################################
/**Reads an unpartitioned mesh from a wavefront .obj file.*/
void chi_mesh::UnpartitionedMesh::ReadFromMsh(const Options &options)
{
  const std::string fname = __FUNCTION__;
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

  //===================================================== Reading every line and
  //                                                      determining size
  std::string file_line;
  std::istringstream iss;
  const std::string node_section_name="$Nodes";
  const std::string elements_section_name = "$Elements";
  const std::string format_section_name = "$MeshFormat";

  //=================================================== Check the format of this input
  while (std::getline(file, file_line))
    if ( format_section_name == file_line ) break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  double format;
  if ( !(iss >> format) )
    throw std::logic_error(fname + ": Failed to read the file format.");
  else if (format != 2.2)
    throw std::logic_error(fname + ": Currently, only msh format 2.2 is supported.");

  //=================================================== Find section with node information
  //                                                    and then read information
  file.seekg(0);
  while (std::getline(file, file_line))
    if ( node_section_name == file_line ) break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  int num_nodes;
  if ( !(iss >> num_nodes) )
    throw std::logic_error(fname + ": Failed while trying to read "
                                   "the number of nodes.");

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
      throw std::logic_error(fname + ": Failed to read vertex index.");

    if (!(iss >> vertices[vert_index-1]->x >> vertices[vert_index-1]->y >> vertices[vert_index-1]->z))
      throw std::logic_error(fname + ": Failed while reading the vertex coordinates.");
  }

  //================================================== Define utility lambdas
  /**Lambda for reading nodes.*/
  auto ReadNodes = [&iss,&fname](int num_nodes)
  {
    std::vector<int> raw_nodes(num_nodes,0);
    for (int i=0; i<num_nodes; ++i)
      if ( !(iss >> raw_nodes[i]) )
        throw std::logic_error(fname + ": Failed when reading element node index.");

    std::vector<uint64_t> nodes(num_nodes,0);
    for (int i=0; i<num_nodes; ++i)
      if ((raw_nodes[i]-1)>=0)
        nodes[i] = raw_nodes[i]-1;
    return nodes;
  };

  /**Lamda for checking if an element is 2D.*/
  auto IsElementType2D = [](int element_type)
  {
    if (element_type == 2 or element_type == 3)
      return true;

    return false;
  };

  /**Lamda for checking if an element is 2D.*/
  auto IsElementType3D = [](int element_type)
  {
    if (element_type >= 4 and element_type <= 7)
      return true;

    return false;
  };

  /**Lambda for checking supported elements.*/
  auto IsElementSupported = [](int element_type)
  {
    if (element_type >= 2 and element_type <= 7)
      return true;

    return false;
  };

  //================================================== Determine mesh type 2D/3D
  bool mesh_is_2D_assumption = true;
  file.seekg(0);
  while (std::getline(file, file_line))
    if ( elements_section_name == file_line ) break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  int num_elems;
  if (!(iss >> num_elems))
    throw std::logic_error(fname + ": Failed to read number of elements.");

  for (int n=0; n<num_elems; n++)
  {
    int elem_type, num_tags, physical_reg, tag, element_index;

    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    if ( !(iss >> element_index >> elem_type >> num_tags) )
      throw std::logic_error(fname + ": Failed while reading element index, "
                                     "element type, and number of tags.");

    if( !(iss>>physical_reg) )
      throw std::logic_error(fname + ": Failed while reading physical region.");

    for (int i=1; i<num_tags; i++)
      if( !(iss >> tag) )
        throw std::logic_error(fname + ": Failed when reading tags.");

    if (IsElementType3D(elem_type))
    {
      mesh_is_2D_assumption = false;
      break; //have the answer now leave loop
    }

    if (elem_type == 15 or elem_type == 1) //skip point and edge
      continue;

    if (not IsElementSupported(elem_type))
      throw std::logic_error(fname + ": Unsupported element encountered.");
  }//for n

  //================================================== Return to the element
  //                                                   listing section
  file.seekg(0);
  while (std::getline(file, file_line))
    if ( elements_section_name == file_line ) break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  if (!(iss >> num_elems))
    throw std::logic_error(fname + ": Failed to read number of elements.");

  for (int n=0; n<num_elems; n++)
  {
    int elem_type, num_tags, physical_reg, tag, element_index;

    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    if ( !(iss >> element_index >> elem_type >> num_tags) )
      throw std::logic_error(fname + ": Failed while reading element index, "
                                     "element type, and number of tags.");

    if( !(iss>>physical_reg) )
      throw std::logic_error(fname + ": Failed while reading physical region.");

    for (int i=1; i<num_tags; i++)
      if ( !(iss >> tag) )
        throw std::logic_error(fname + ": Failed when reading tags.");

    if (elem_type == 15 or elem_type == 1) //skip point and edge type elements
      continue;

    if (not mesh_is_2D_assumption and IsElementType2D(elem_type))
      continue;

    if (mesh_is_2D_assumption and (not IsElementType2D(elem_type)))
      continue;

    if (not IsElementSupported(elem_type))
      throw std::logic_error(fname + ": Unsupported element encountered.");

    int num_cell_nodes;
    if (elem_type == 2)      //3-node triangle
      num_cell_nodes = 3;
    else if (elem_type == 3 or elem_type == 4) //4-node quadrangle or tet
      num_cell_nodes = 4;
    else if (elem_type == 5) //8-node hexahedron
      num_cell_nodes = 8;
    else
      continue;

    //====================================== Make the cell
    if (mesh_is_2D_assumption)
      raw_cells.push_back(new LightWeightCell(CellType::POLYGON));
    else
      raw_cells.push_back(new LightWeightCell(CellType::POLYHEDRON));

    auto& cell = *raw_cells.back();
    cell.material_id = physical_reg;
    cell.vertex_ids = ReadNodes(num_cell_nodes);

    //====================================== Populate faces
    if (elem_type == 2 or elem_type == 3) //3-node triangle or 4-node quadrangle
    {
      size_t num_verts = cell.vertex_ids.size();
      for (int e=0; e<num_verts; e++)
      {
        int ep1 = (e<(num_verts-1))? e+1 : 0;
        LightWeightFace face;

        face.vertex_ids = {cell.vertex_ids[e], cell.vertex_ids[ep1]};

        cell.faces.push_back(std::move(face));
      }//for e
    }//if 2D elements
    else if (elem_type == 4) // 4-node tetrahedron
    {
      auto& v = cell.vertex_ids;
      std::vector<LightWeightFace> lw_faces(4);
      lw_faces[0].vertex_ids = {v[0], v[2], v[1]};
      lw_faces[1].vertex_ids = {v[0], v[3], v[2]};
      lw_faces[2].vertex_ids = {v[3], v[1], v[2]};
      lw_faces[3].vertex_ids = {v[3], v[0], v[1]};

      for (auto& lw_face : lw_faces) cell.faces.push_back(lw_face);
    }
    else if (elem_type == 5) //8-node hexahedron
    {
      auto& v = cell.vertex_ids;
      std::vector<LightWeightFace> lw_faces(6);
      lw_faces[0].vertex_ids = {v[5], v[1], v[2], v[6]}; //East face
      lw_faces[1].vertex_ids = {v[0], v[4], v[7], v[3]}; //West face
      lw_faces[2].vertex_ids = {v[0], v[3], v[2], v[1]}; //North face
      lw_faces[3].vertex_ids = {v[4], v[5], v[6], v[7]}; //South face
      lw_faces[4].vertex_ids = {v[2], v[3], v[7], v[6]}; //Top face
      lw_faces[5].vertex_ids = {v[0], v[1], v[5], v[4]}; //Bottom face

      for (auto& lw_face : lw_faces) cell.faces.push_back(lw_face);
    }

  }//for elements

  file.close();

  //======================================== Always do this
  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();

  chi_log.Log() << "Done processing " << options.file_name << ".\n"
                << "Number of nodes read: " << vertices.size() << "\n"
                << "Number of cells read: " << raw_cells.size();
}



