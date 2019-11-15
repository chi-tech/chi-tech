#include"chi_surfacemesh.h"
#include<iostream>
#include<fstream>
#include <algorithm>

#include <chi_log.h>

extern ChiLog chi_log;

//#########################################################
/** Loads a surface mesh from a wavefront .obj file.*/
int chi_mesh::SurfaceMesh::
    ImportFromOBJFile(const char* fileName, bool as_poly=false)
{


  //===================================================== Opening the file
  std::ifstream file;
  file.open(fileName);
  if (!file.is_open())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Failed to open file: "<< fileName<<" in call "
      << "to ImportFromOBJFile \n";
    exit(EXIT_FAILURE);
  }

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
    if (first_word.compare("v")==0) {
      chi_mesh::Vertex* newVertex = new chi_mesh::Vertex;
      for (int k=1;k<=3;k++)
      {
        //================================== Extract sub word
        beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
        end_of_word = file_line.find(delimiter, beg_of_word);
        sub_word = file_line.substr(beg_of_word, end_of_word-beg_of_word);

        //============================= Convert word to number
        try{ double numValue = std::stod(sub_word);
          //*newVertex->value[k-1] = numValue;
          if (k==1)
          {
            newVertex->x = numValue;
          }
          else if (k==2)
          {
            newVertex->y = numValue;
          }
          else if (k==3)
          {
            newVertex->z = numValue;
          }
        }

        //============================= Catch conversion error
        catch(const std::invalid_argument& ia)
        {
          //*newVertex->value[k-1] = 0.0;
          std::cout<<"Exception caught!"<<std::endl;
        }

        //============================= Stop word extraction on line end
        if (end_of_word==std::string::npos) {break;}
      }
      this->vertices.push_back(*newVertex);
    }

    //===================================================== Keyword "vt" for Vertex
    if (first_word.compare("vt")==0) {
      chi_mesh::Vertex* newVertex = new chi_mesh::Vertex;
      for (int k=1;k<=2;k++)
      {
        //================================== Extract sub word
        beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
        end_of_word = file_line.find(delimiter, beg_of_word);
        sub_word = file_line.substr(beg_of_word, end_of_word-beg_of_word);

        //============================= Convert word to number
        try{ double numValue = std::stod(sub_word);
          //*newVertex->value[k-1] = numValue;
          if (k==1)
          {
            newVertex->x = numValue;
          }
          else if (k==2)
          {
            newVertex->y = numValue;
          }
          else if (k==3)
          {
            newVertex->z = numValue;
          }
        }

        //============================= Catch conversion error
        catch(const std::invalid_argument& ia)
        {
          //*newVertex->value[k-1] = 0.0;
          std::cout<<"Exception caught!"<<std::endl;
        }

        //============================= Stop word extraction on line end
        if (end_of_word==std::string::npos) {break;}
      }
      this->tex_vertices.push_back(*newVertex);
    }

    //===================================================== Keyword "vn" for normal
    if (first_word.compare("vn")==0) {
      chi_mesh::Normal* newNormal = new chi_mesh::Normal;
      for (int k=1;k<=3;k++)
      {
        //================================== Extract sub word
        beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
        end_of_word = file_line.find(delimiter, beg_of_word);
        sub_word = file_line.substr(beg_of_word, end_of_word-beg_of_word);

        //============================= Convert word to number
        try{ double numValue = std::stod(sub_word);
          //*newNormal->value[k-1] = numValue;
          if (k==1)
          {
            newNormal->x = numValue;
          }
          else if (k==2)
          {
            newNormal->y = numValue;
          }
          else if (k==3)
          {
            newNormal->z = numValue;
          }
        }

        //============================= Catch conversion error
        catch(const std::invalid_argument& ia)
        {
          //*newNormal->value[k-1] = 0.0;
          std::cout<<"Exception caught!"<<std::endl;
        }

        //============================= Stop word extraction on line end
        if (end_of_word==std::string::npos) {break;}
      }
      this->normals.push_back(*newNormal);
    }

    //===================================================== Keyword "f" for face
    if (first_word.compare("f")==0)
    {
      int number_of_verts = std::count(file_line.begin(),file_line.end(),'/')/2;
      if ((number_of_verts==3) && (!as_poly))
      {
        chi_mesh::Face* newFace = new chi_mesh::Face;

        for (int k=1;k<=3;k++)
        {
          //================================== Extract sub word
          beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
          end_of_word = file_line.find(delimiter, beg_of_word);
          sub_word = file_line.substr(beg_of_word, end_of_word-beg_of_word);

          //============================= Extract locations of hiphens
          size_t first_dash = sub_word.find("/");
          size_t last_dash  = sub_word.find_last_of("/");

          //============================= Extract the words ass. w vertex and normal
          std::string vert_word = sub_word.substr(0,first_dash-0);
          std::string norm_word = sub_word.substr(last_dash+1,sub_word.length()-last_dash-1);

          //============================= Convert word to number (Vertex)
          try{ int numValue = std::stoi(vert_word);
            newFace->v_index[k-1] = numValue-1;
          }
          catch(const std::invalid_argument& ia){std::cout<<"Exception caught!"<<std::endl; }

          //============================= Convert word to number (Normal)
          try{ int numValue = std::stoi(norm_word);
            newFace->n_index[k-1] = numValue-1;
          }
          catch(const std::invalid_argument& ia){std::cout<<"Exception caught!"<<std::endl; }

          //============================= Convert word to number (Texture Vertex)
          if (last_dash>(first_dash+1)){
            std::string tvert_word = sub_word.substr(first_dash+1,last_dash-first_dash-1);
            try{ int numValue = std::stoi(tvert_word);
              newFace->vt_index[k-1] = numValue-1;
            }
            catch(const std::invalid_argument& ia)
            {
              std::cout<<"Exception caught!"<<std::endl;
            }
          }

          //============================= Stop word extraction on line end
          if (end_of_word==std::string::npos) {break;}
        }


        //==================================== Set edges
        newFace->e_index[0][0] = newFace->v_index[0];
        newFace->e_index[0][1] = newFace->v_index[1];

        newFace->e_index[1][0] = newFace->v_index[1];
        newFace->e_index[1][1] = newFace->v_index[2];

        newFace->e_index[2][0] = newFace->v_index[2];
        newFace->e_index[2][1] = newFace->v_index[0];

        this->faces.push_back(*newFace);
      }
      else
      {
        chi_mesh::PolyFace* newFace = new chi_mesh::PolyFace;

        for (int k=1;k<=number_of_verts;k++)
        {
          //================================== Extract sub word
          beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
          end_of_word = file_line.find(delimiter, beg_of_word);
          sub_word = file_line.substr(beg_of_word, end_of_word-beg_of_word);

          //============================= Extract locations of hiphens
          size_t first_dash = sub_word.find("/");
          size_t last_dash  = sub_word.find_last_of("/");

          //============================= Extract the words ass. w vertex and normal
          std::string vert_word = sub_word.substr(0,first_dash-0);
          std::string norm_word = sub_word.substr(last_dash+1,sub_word.length()-last_dash-1);

          //============================= Convert word to number (Vertex)
          try{ int numValue = std::stoi(vert_word);
            newFace->v_indices.push_back(numValue-1);
          }
          catch(const std::invalid_argument& ia){std::cout<<"Exception caught!"<<std::endl; }

          //============================= Convert word to number (Normal)
          try{ int numValue = std::stoi(norm_word);
            //newFace->n_indices.push_back(numValue-1);
          }
          catch(const std::invalid_argument& ia){std::cout<<"Exception caught!"<<std::endl; }

          //============================= Convert word to number (Texture Vertex)
          if (last_dash>(first_dash+1)){
            std::string tvert_word = sub_word.substr(first_dash+1,last_dash-first_dash-1);
            try{ int numValue = std::stoi(tvert_word);
              //newFace->vt_indices.push_back(numValue-1);
            }
            catch(const std::invalid_argument& ia)
            {
              std::cout<<"Exception caught!"<<std::endl;
            }
          }

          //============================= Stop word extraction on line end
          if (end_of_word==std::string::npos) {break;}
        }

        for (int v=0;v<(newFace->v_indices.size());v++)
        {
          int* side_indices = new int[4];

          side_indices[0] = newFace->v_indices[v];
          side_indices[1] = newFace->v_indices[v+1];
          side_indices[2] = -1;
          side_indices[3] = -1;

          if ((v+1)>=newFace->v_indices.size())
          {
            side_indices[1] = newFace->v_indices[0];
          }

          newFace->edges.push_back(side_indices);
        }


        this->poly_faces.push_back(newFace);
      }


    }


    //=================================================== Keyword "l" for line
    if (first_word.compare("l")==0){
      chi_mesh::Edge newEdge;

      //================================== Extract sub word
      beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
      end_of_word = file_line.find(delimiter, beg_of_word);
      sub_word = file_line.substr(beg_of_word, end_of_word-beg_of_word);

      //================================== Convert to number
      try{ int numValue = std::stoi(sub_word);
        newEdge.v_index[0] = numValue-1;
      }
      catch(const std::invalid_argument& ia){std::cout<<"Exception caught!"<<std::endl; }

      //================================== Extract sub word
      beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
      end_of_word = file_line.find(delimiter, beg_of_word);
      sub_word = file_line.substr(beg_of_word, end_of_word-beg_of_word);

      //================================== Convert to number
      try{ int numValue = std::stoi(sub_word);
        newEdge.v_index[1] = numValue-1;
      }
      catch(const std::invalid_argument& ia){std::cout<<"Exception caught!"<<std::endl; }

      newEdge.vertices[0] = this->vertices.at(newEdge.v_index[0]);
      newEdge.vertices[1] = this->vertices.at(newEdge.v_index[1]);

      this->lines.push_back(newEdge);

      //printf("line %d->%d\n",newEdge.v_index[0],newEdge.v_index[1]);
    }

  }
  file.close();

  //======================================================= Calculate face properties
  std::vector<chi_mesh::Face>::iterator curFace;
  for (curFace = this->faces.begin(); curFace!=this->faces.end(); curFace++)
  {
    //=========================================== Calculate geometrical normal
    chi_mesh::Vertex vA = this->vertices.at(curFace->v_index[0]);
    chi_mesh::Vertex vB = this->vertices.at(curFace->v_index[1]);
    chi_mesh::Vertex vC = this->vertices.at(curFace->v_index[2]);

    chi_mesh::Vector vAB = vB-vA;
    chi_mesh::Vector vBC = vC-vB;

    curFace->geometric_normal = vAB.Cross(vBC);
    curFace->geometric_normal = curFace->geometric_normal/curFace->geometric_normal.Norm();

    //=========================================== Calculate Assigned normal
    chi_mesh::Vertex nA = this->normals.at(curFace->n_index[0]);
    chi_mesh::Vertex nB = this->normals.at(curFace->n_index[1]);
    chi_mesh::Vertex nC = this->normals.at(curFace->n_index[2]);

    chi_mesh::Vector nAvg = (nA+nB+nC)/3.0;
    nAvg = nAvg/nAvg.Norm();

    curFace->assigned_normal = nAvg;

    //=========================================== Compute face center
    curFace->face_centroid = (vA+vB+vC)/3.0;
  }
  std::vector<chi_mesh::PolyFace*>::iterator curPFace;
  for (curPFace = this->poly_faces.begin();
       curPFace!=this->poly_faces.end();
       curPFace++)
  {
    chi_mesh::Vector centroid;
    int num_verts = (*curPFace)->v_indices.size();
    for (int v=0; v<num_verts; v++)
      centroid = centroid + vertices[(*curPFace)->v_indices[v]];

    centroid = centroid/num_verts;

    (*curPFace)->face_centroid = centroid;

    chi_mesh::Vector n = (vertices[(*curPFace)->v_indices[1]] -
                          vertices[(*curPFace)->v_indices[0]]).Cross(
                          centroid - vertices[(*curPFace)->v_indices[1]]);
    n = n/n.Norm();

    (*curPFace)->geometric_normal = n;
  }

  UpdateInternalConnectivity();

  //============================================= Check each vertex is accounted
  chi_log.Log(LOG_0)
  << "Surface mesh loaded with "
  << this->faces.size() << " triangle faces and "
  << this->poly_faces.size() << " polygon faces.";
  //exit(EXIT_FAILURE);

  return 0;
}

//#########################################################
/** Loads a surface mesh from triangle's file format.*/
int chi_mesh::SurfaceMesh::
ImportFromTriangleFiles(const char* fileName, bool as_poly=false)
{
  std::string node_filename = std::string(fileName) +
                              std::string(".1.node");
  std::string tria_filename = std::string(fileName) +
                              std::string(".1.ele");

  //===================================================== Opening the node file
  std::ifstream file;
  file.open(node_filename);
  if (!file.is_open())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Failed to open file: "<< node_filename <<" in call "
      << "to ImportFromOBJFile \n";
    exit(EXIT_FAILURE);
  }

  int num_verts;
  char line[250];
  file >> num_verts;
  file.getline(line,250);
  for (int v=1; v<=num_verts; v++)
  {
    int vert_index;
    chi_mesh::Vertex vertex;
    file >> vert_index >> vertex.x >> vertex.y;
    file.getline(line,250);

    vertices.push_back(vertex);
  }

  file.close();

  //===================================================== Opening the ele file
  file.open(tria_filename);
  if (!file.is_open())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Failed to open file: "<< tria_filename <<" in call "
      << "to ImportFromOBJFile \n";
    exit(EXIT_FAILURE);
  }

  int num_tris;

  file >> num_tris;
  file.getline(line,250);
  for (int v=1; v<=num_tris; v++)
  {
    int tri_index;
    chi_mesh::PolyFace* newFace = new chi_mesh::PolyFace;

    int v0,v1,v2;
    file >> tri_index >> v0 >> v1 >> v2;
    file.getline(line,250);


    newFace->v_indices.resize(3);
    newFace->v_indices[0] = v0-1;
    newFace->v_indices[1] = v1-1;
    newFace->v_indices[2] = v2-1;

    for (int e=0; e<3; e++)
    {
      int* side_indices = new int[4];

      if (e<2)
      {
        side_indices[0] = newFace->v_indices[e];
        side_indices[1] = newFace->v_indices[e+1];
        side_indices[2] = -1;
        side_indices[3] = -1;
      }
      else
      {
        side_indices[0] = newFace->v_indices[e];
        side_indices[1] = newFace->v_indices[0];
        side_indices[2] = -1;
        side_indices[3] = -1;
      }
      newFace->edges.push_back(side_indices);
    }

    poly_faces.push_back(newFace);
  }

  file.close();

  //======================================================= Calculate face properties
  std::vector<chi_mesh::Face>::iterator curFace;
  for (curFace = this->faces.begin(); curFace!=this->faces.end(); curFace++)
  {
    //=========================================== Calculate geometrical normal
    chi_mesh::Vertex vA = this->vertices.at(curFace->v_index[0]);
    chi_mesh::Vertex vB = this->vertices.at(curFace->v_index[1]);
    chi_mesh::Vertex vC = this->vertices.at(curFace->v_index[2]);

    chi_mesh::Vector vAB = vB-vA;
    chi_mesh::Vector vBC = vC-vB;

    curFace->geometric_normal = vAB.Cross(vBC);
    curFace->geometric_normal = curFace->geometric_normal/curFace->geometric_normal.Norm();

    //=========================================== Calculate Assigned normal
    chi_mesh::Vertex nA = this->normals.at(curFace->n_index[0]);
    chi_mesh::Vertex nB = this->normals.at(curFace->n_index[1]);
    chi_mesh::Vertex nC = this->normals.at(curFace->n_index[2]);

    chi_mesh::Vector nAvg = (nA+nB+nC)/3.0;
    nAvg = nAvg/nAvg.Norm();

    curFace->assigned_normal = nAvg;

    //=========================================== Compute face center
    curFace->face_centroid = (vA+vB+vC)/3.0;
  }
  std::vector<chi_mesh::PolyFace*>::iterator curPFace;
  for (curPFace = this->poly_faces.begin();
       curPFace!=this->poly_faces.end();
       curPFace++)
  {
    chi_mesh::Vector centroid;
    int num_verts = (*curPFace)->v_indices.size();
    for (int v=0; v<num_verts; v++)
      centroid = centroid + vertices[(*curPFace)->v_indices[v]];

    centroid = centroid/num_verts;

    (*curPFace)->face_centroid = centroid;

    chi_mesh::Vector n = (vertices[(*curPFace)->v_indices[1]] -
                          vertices[(*curPFace)->v_indices[0]]).Cross(
      centroid - vertices[(*curPFace)->v_indices[1]]);
    n = n/n.Norm();

    (*curPFace)->geometric_normal = n;
  }

  UpdateInternalConnectivity();

  //============================================= Check each vertex is accounted
  chi_log.Log(LOG_0)
    << "Surface mesh loaded with "
    << this->faces.size() << " triangle faces and "
    << this->poly_faces.size() << " polygon faces.";
  //exit(EXIT_FAILURE);

  return 0;
}



//#########################################################
/** Runs over the faces of the surfacemesh and determines
 * neighbors. */
void chi_mesh::SurfaceMesh::UpdateInternalConnectivity()
{
  std::vector<chi_mesh::Face>::iterator curFace;
  for (curFace = this->faces.begin(); curFace!=this->faces.end(); curFace++)
  {
    int outer_index = std::distance(this->faces.begin(),curFace);

    //Reset face connections
    for (int e=0;e<3;e++) //Loop over current face's edges
    {
      curFace->e_index[e][2] = -1; //tri index
      curFace->e_index[e][3] = -1;
    }
    //printf("Outer index %d\n",outer_index);
    //=========================================== Determine face connectivity
    //If the head v_index of the current triangle is equal to the tail
    //of the other triangle and vice versa on the tail of the current
    //triangle, then they are connected.
    std::vector<chi_mesh::Face>::iterator other_face;
    for (other_face = faces.begin(); other_face!=faces.end(); other_face++)
    {
      int inner_index = std::distance(this->faces.begin(),other_face);
      //if (outer_index!=inner_index)
      {
        for (int e=0;e<3;e++) //Loop over current face's edges
        {

          for (int e2=0;e2<3;e2++)  //Loop over other face's edges
          {
            if ( (curFace->e_index[e][0]==other_face->e_index[e2][1]) &&
                 (curFace->e_index[e][1]==other_face->e_index[e2][0]) )
            {
              curFace->e_index[e][2] = inner_index; //tri index
              curFace->e_index[e][3] = e2;                       //edge index
            }
          }
        }
      }

    }
  }

  for (int f=0; f<this->poly_faces.size(); f++)
  {
    chi_mesh::PolyFace* curFace = this->poly_faces[f];

    for (int f2=0; f2<this->poly_faces.size(); f2++)
    {
      if (f!=f2)
      {
        chi_mesh::PolyFace* other_face = this->poly_faces[f2];
        for (int e=0;e<curFace->edges.size();e++) //Loop over current face's edges
        {

          for (int e2=0;e2<other_face->edges.size();e2++)  //Loop over other face's edges
          {
            if ( (curFace->edges[e][0]==other_face->edges[e2][1]) &&
                 (curFace->edges[e][1]==other_face->edges[e2][0]) )
            {
              curFace->edges[e][2] = f2; //tri index
              curFace->edges[e][3] = e2;                       //edge index
            }
          }
        }
      }
    }
  }
}



//#########################################################
/**Exports the triangular faces of a surface mesh to
 * wavefront .obj files.*/
void chi_mesh::SurfaceMesh::ExportToOBJFile(const char *fileName)
{

//  if (this->faces.empty())
//  {
//    std::cout << "Cannot export empty SurfaceMesh\n";
//    return;
//  }
  FILE* outputFile = fopen(fileName,"w");
  if (outputFile==NULL)
  {
    printf("Error creating file %s!\n",fileName);
    return;
  }

  fprintf(outputFile,"# Exported mesh file from tringulation script\n");
  fprintf(outputFile,"o %s\n","ChitechTriMesh");

  std::vector<chi_mesh::Vertex>::iterator cur_v;
  for (cur_v = this->vertices.begin();
       cur_v != this->vertices.end();
       cur_v++)
  {
    fprintf(outputFile,"v %9.6f %9.6f %9.6f\n",cur_v->x,cur_v->y,cur_v->z);
  }

  for (unsigned ell=0; ell<this->lines.size(); ell++)
  {
    fprintf(outputFile,"l %d %d \n",lines[ell].v_index[0]+1,
                                    lines[ell].v_index[1]+1);
  }

  if (faces.size()>0)
  {
    chi_mesh::Face first_face = this->faces.front();
    fprintf(outputFile,"vn %.4f %.4f %.4f\n", first_face.geometric_normal.x,
            first_face.geometric_normal.y,
            first_face.geometric_normal.z);
    fprintf(outputFile,"s off\n");

    std::vector<chi_mesh::Face>::iterator cur_face;
    for (cur_face = this->faces.begin();
         cur_face != this->faces.end();
         cur_face++)
    {
      fprintf(outputFile,"f %d//1 %d//1 %d//1\n",cur_face->v_index[0]+1,
              cur_face->v_index[1]+1,
              cur_face->v_index[2]+1);
    }
  }
  if (poly_faces.size()>0)
  {
    chi_mesh::PolyFace* first_face = this->poly_faces.front();
    fprintf(outputFile,"vn %.4f %.4f %.4f\n", first_face->geometric_normal.x,
            first_face->geometric_normal.y,
            first_face->geometric_normal.z);
    fprintf(outputFile,"s off\n");

    for (int f=0; f<poly_faces.size(); f++)
    {
      fprintf(outputFile,"f ");
      for (int v=0; v<poly_faces[f]->v_indices.size(); v++)
      {
        fprintf(outputFile,"%d//1 ",poly_faces[f]->v_indices[v]+1);
      }
      fprintf(outputFile,"\n");
    }
  }


  fclose(outputFile);
  printf("Exported mesh to %s\n",fileName);
}

//#########################################################
/**Exports a PSLG to triangle1.6's .poly format.*/
void chi_mesh::SurfaceMesh::ExportToPolyFile(const char *fileName)
{
  FILE* outputFile = fopen(fileName,"w");
  if (outputFile==NULL)
  {
    printf("Error creating file %s!\n",fileName);
    return;
  }

  fprintf(outputFile,"%lu 2 0 0\n", vertices.size());
  for (int v=0; v<vertices.size(); v++)
  {
    fprintf(outputFile,"%d %.15f %.15f 0\n",v+1,vertices[v].x,vertices[v].y);
  }

  fprintf(outputFile,"%lu 0\n", lines.size());
  for (int e=0; e<lines.size(); e++)
  {
    fprintf(outputFile,"%d %d %d\n",e+1,lines[e].v_index[0]+1,
                                      lines[e].v_index[1]+1);
  }

  fprintf(outputFile,"0");


  fclose(outputFile);
  printf("Exported mesh to %s\n",fileName);
}
