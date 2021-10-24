#include "chi_unpartitioned_mesh.h"

#include "chi_log.h"
extern ChiLog& chi_log;
#include "chi_mpi.h"

#include <algorithm>

//###################################################################
/**Reads an unpartitioned mesh from a wavefront .obj file.*/
void chi_mesh::UnpartitionedMesh::ReadFromWavefrontOBJ(const Options &options)
{
  //===================================================== Opening the file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Failed to open file: "<< options.file_name<<" in call "
      << "to ImportFromOBJFile \n";
    exit(EXIT_FAILURE);
  }

  chi_log.Log() << "Making Unpartitioned mesh from wavefront file "
                << options.file_name;
  MPI_Barrier(MPI_COMM_WORLD);

  //===================================================== Reading every line and determining size
  std::string file_line;
  std::string delimiter = " ";
  int counter=0;
  while (std::getline(file,file_line))
  {
    counter++;

    //===================================================== Get the first word
    size_t beg_of_word = file_line.find_first_not_of(delimiter);
    size_t end_of_word = file_line.find(delimiter,beg_of_word-beg_of_word);
    std::string first_word = file_line.substr(beg_of_word,end_of_word);
    std::string sub_word;

    //===================================================== Keyword "v" for Vertex
    if (first_word == "v")
    {
      chi_mesh::Vertex newVertex;
      for (int k=1;k<=3;k++)
      {
        //================================== Extract sub word
        beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
        end_of_word = file_line.find(delimiter, beg_of_word);
        sub_word = file_line.substr(beg_of_word, end_of_word-beg_of_word);

        //============================= Convert word to number
        try
        {
          double numValue = std::stod(sub_word);

          if (k==1)
            newVertex.x = numValue;
          else if (k==2)
            newVertex.y = numValue;
          else if (k==3)
            newVertex.z = numValue;
        }

          //============================= Catch conversion error
        catch(const std::invalid_argument& ia)
        {
          std::cout<<"Exception caught!"<<std::endl;
        }

        //============================= Stop word extraction on line end
        if (end_of_word==std::string::npos) {break;}
      }
      this->vertices.push_back(newVertex);
    }

    //===================================================== Keyword "f" for face
    if (first_word == "f")
    {
      size_t number_of_verts =
        std::count(file_line.begin(), file_line.end(), '/')/2;

      CellType sub_type = CellType::POLYGON;
      if      (number_of_verts == 3) sub_type = CellType::TRIANGLE;
      else if (number_of_verts == 4) sub_type = CellType::QUADRILATERAL;

      auto cell = new LightWeightCell(CellType::POLYGON, sub_type);

      for (size_t k=1;k<=number_of_verts;k++)
      {
        //================================== Extract sub word
        beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
        end_of_word = file_line.find(delimiter, beg_of_word);
        sub_word = file_line.substr(beg_of_word, end_of_word-beg_of_word);

        //============================= Extract locations of hiphens
        size_t first_dash = sub_word.find('/');
        size_t last_dash  = sub_word.find_last_of('/');

        //============================= Extract the words ass. w vertex and normal
        std::string vert_word = sub_word.substr(0,first_dash-0);
        std::string norm_word = sub_word.substr(last_dash+1,sub_word.length()-last_dash-1);

        //============================= Convert word to number (Vertex)
        try{ int numValue = std::stoi(vert_word);
          cell->vertex_ids.push_back(numValue-1);
        }
        catch(const std::invalid_argument& ia){std::cout<<"Exception caught!"<<std::endl; }

        //============================= Stop word extraction on line end
        if (end_of_word==std::string::npos) {break;}
      }

      size_t num_verts = cell->vertex_ids.size();
      for (uint64_t v=0;v<num_verts; ++v)
      {
        LightWeightFace face;

        face.vertex_ids.resize(2);
        face.vertex_ids[0] = cell->vertex_ids[v];
        face.vertex_ids[1] = (v<(num_verts-1))? cell->vertex_ids[v+1] :
                                                cell->vertex_ids[0];

        cell->faces.push_back(std::move(face));
      }

      this->raw_cells.push_back(cell);
    }
  }
  file.close();

  //======================================== Always do this
  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();
}