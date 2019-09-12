#include "chi_meshcontinuum.h"
#include <fstream>
#include "../Cell/cell_polyhedron.h"
#include "../Cell/cell_polygon.h"
#include <ChiPhysics/chi_physics.h>
#include <ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h>


extern ChiPhysics chi_physics_handler;

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

//###################################################################
/**Export cells to python.*/
void chi_mesh::MeshContinuum::
ExportCellsToObj(const char* fileName, bool per_material,
                    int options)
{
  if (!per_material)
  {
    FILE* of = fopen(fileName,"w");

    if (of==NULL)
    {
      chi_log.Log(LOG_ALLERROR) << "Could not open file: "
                                  << std::string(fileName);
      exit(EXIT_FAILURE);
    }

    //====================================== Develop list of faces and nodes
    std::set<int> nodes_set;
    std::vector<chi_mesh::PolyFace*> faces_to_export;
    for (int c=0; c<local_cell_glob_indices.size(); c++)
    {
      int cell_glob_index = local_cell_glob_indices[c];

      if (dynamic_cast<chi_mesh::CellPolyhedron*>(cells[c]))
      {
        chi_mesh::CellPolyhedron* polyh_cell =
          ((chi_mesh::CellPolyhedron*)cells[c]);

        for (int f=0; f<polyh_cell->faces.size(); f++)
        {
          if (polyh_cell->faces[f]->face_indices[NEIGHBOR]<0)
          {
            faces_to_export.push_back(polyh_cell->faces[f]);

            for (int v=0; v<polyh_cell->faces[f]->v_indices.size(); v++)
            {
              nodes_set.insert(polyh_cell->faces[f]->v_indices[v]);
            }
          }//if boundary
        }//for face
      }//if polyhedron
    }//for local cell

    //====================================== Write header
    fprintf(of,"# Exported mesh file from Extrusion script\n");
    std::string str_file_name(fileName);
    std::string file_base_name =
      str_file_name.substr(0,str_file_name.find("."));
    fprintf(of,"o %s\n",file_base_name.c_str());

    //====================================== Develop node mapping and write them
    std::vector<int> node_mapping(nodes.size(),-1);
    std::set<int>::iterator node;
    int node_counter=0;
    for (node =  nodes_set.begin();
         node != nodes_set.end();
         node++)
    {
      node_counter++;
      int node_g_index = *node;
      node_mapping[node_g_index] = node_counter;

      chi_mesh::Vertex* cur_v = nodes[node_g_index];

      fprintf(of,"v %9.6f %9.6f %9.6f\n",cur_v->x,cur_v->y,cur_v->z);
    }

    //====================================== Write face normals
    for (int f=0; f<faces_to_export.size(); f++)
    {
      fprintf(of,"vn %.4f %.4f %.4f\n",
              faces_to_export[f]->geometric_normal.x,
              faces_to_export[f]->geometric_normal.y,
              faces_to_export[f]->geometric_normal.z);
    }

    //====================================== Write faces
    int normal_counter=0;
    for (int f=0; f<faces_to_export.size(); f++)
    {
      normal_counter++;
      fprintf(of,"f");
      for (int v=0; v<faces_to_export[f]->v_indices.size(); v++)
      {
        int v_g_index = faces_to_export[f]->v_indices[v];
        int v_mapped  = node_mapping[v_g_index];

        fprintf(of," %d//%d",v_mapped,normal_counter);
      }
      fprintf(of,"\n");
    }


    fclose(of);

    chi_log.Log(LOG_0)
     << "Exported Volume mesh to "
     << str_file_name;
  }//Whole mesh
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PER MATERIAL
  else
  {
    //========================================= Get base name
    std::string str_file_name(fileName);
    std::string file_base_name =
      str_file_name.substr(0,str_file_name.find("."));

    if (chi_physics_handler.material_stack.size() == 0)
    {
      chi_log.Log(LOG_0WARNING)
        << "ExportCellsToObj: No mesh will be exported because there "
        << "are no physics materials present";
    }


    for (int mat=0; mat<chi_physics_handler.material_stack.size(); mat++)
    {
      std::string mat_base_name = file_base_name +
                                  std::string("_m") +
                                  std::to_string(mat);
      std::string mat_file_name = mat_base_name +
                                  std::string(".obj");
      FILE* of = fopen(mat_file_name.c_str(),"w");

      if (of==NULL)
      {
        chi_log.Log(LOG_ALLERROR) << "Could not open file: "
                                  << mat_file_name;
        exit(EXIT_FAILURE);
      }

      //====================================== Develop list of faces and nodes
      std::set<int> nodes_set;
      std::vector<chi_mesh::PolyFace*> faces_to_export;
      for (int c=0; c<local_cell_glob_indices.size(); c++)
      {
        int cell_glob_index = local_cell_glob_indices[c];

        if (dynamic_cast<chi_mesh::CellPolyhedron*>(cells[c]))
        {
          chi_mesh::CellPolyhedron* polyh_cell =
            ((chi_mesh::CellPolyhedron*)cells[c]);

          if (polyh_cell->material_id != mat) continue;

          for (int f=0; f<polyh_cell->faces.size(); f++)
          {
            int adjcell_glob_index =
              polyh_cell->faces[f]->face_indices[NEIGHBOR];

            if (adjcell_glob_index<0)
            {
              faces_to_export.push_back(polyh_cell->faces[f]);

              for (int v=0; v<polyh_cell->faces[f]->v_indices.size(); v++)
              {
                nodes_set.insert(polyh_cell->faces[f]->v_indices[v]);
              }
            }//if boundary
            else
            {
              auto adj_cell = cells[adjcell_glob_index];

              if (adj_cell->material_id != mat)
              {
                faces_to_export.push_back(polyh_cell->faces[f]);

                for (int v=0; v<polyh_cell->faces[f]->v_indices.size(); v++)
                {
                  nodes_set.insert(polyh_cell->faces[f]->v_indices[v]);
                }
              }//if material missmatch
            }//if neigbor cell
          }//for face
        }//if polyhedron
      }//for local cell

      //====================================== Write header
      fprintf(of,"# Exported mesh file from Extrusion script\n");
      fprintf(of,"o %s\n",mat_base_name.c_str());

      //====================================== Develop node mapping and write them
      std::vector<int> node_mapping(nodes.size(),-1);
      std::set<int>::iterator node;
      int node_counter=0;
      for (node =  nodes_set.begin();
           node != nodes_set.end();
           node++)
      {
        node_counter++;
        int node_g_index = *node;
        node_mapping[node_g_index] = node_counter;

        chi_mesh::Vertex* cur_v = nodes[node_g_index];

        fprintf(of,"v %9.6f %9.6f %9.6f\n",cur_v->x,cur_v->y,cur_v->z);
      }

      //====================================== Write face normals
      for (int f=0; f<faces_to_export.size(); f++)
      {
        fprintf(of,"vn %.4f %.4f %.4f\n",
                faces_to_export[f]->geometric_normal.x,
                faces_to_export[f]->geometric_normal.y,
                faces_to_export[f]->geometric_normal.z);
      }

      //====================================== Write faces
      int normal_counter=0;
      for (int f=0; f<faces_to_export.size(); f++)
      {
        normal_counter++;
        fprintf(of,"f");
        for (int v=0; v<faces_to_export[f]->v_indices.size(); v++)
        {
          int v_g_index = faces_to_export[f]->v_indices[v];
          int v_mapped  = node_mapping[v_g_index];

          fprintf(of," %d//%d",v_mapped,normal_counter);
        }
        fprintf(of,"\n");
      }


      fclose(of);

      chi_log.Log(LOG_0)
        << "Exported Material Volume mesh to "
        << mat_file_name;
    }//for mat
  }//if per material

}