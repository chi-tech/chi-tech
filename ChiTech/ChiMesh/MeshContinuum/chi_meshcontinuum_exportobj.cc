#include "chi_meshcontinuum.h"
#include <fstream>
#include <ChiPhysics/chi_physics.h>
#include <ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h>


extern ChiPhysics&  chi_physics_handler;

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

//###################################################################
/**Export cells to python.
 *
 * \todo Export Cells to OBJ needs polygon support. */
void chi_mesh::MeshContinuum::
 ExportCellsToObj(const char* fileName, bool per_material, int options) const
{
  if (!per_material)
  {
    FILE* of = fopen(fileName,"w");

    if (of == nullptr)
    {
      chi_log.Log(LOG_ALLERROR) << "Could not open file: "
                                  << std::string(fileName);
      exit(EXIT_FAILURE);
    }

    //====================================== Develop list of faces and nodes
    std::set<int> nodes_set;
    std::vector<chi_mesh::CellFace> faces_to_export;
    for (auto& cell : local_cells)
    {
      if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
      {
        for (auto& face : cell.faces)
        {
          if (not face.has_neighbor)
          {
            faces_to_export.push_back(face);

            for (int vid : face.vertex_ids)
              nodes_set.insert(vid);
          }//if boundary
        }//for face
      }//if polyhedron
    }//for local cell

    //====================================== Write header
    fprintf(of,"# Exported mesh file from Extrusion script\n");
    std::string str_file_name(fileName);
    std::string file_base_name =
      str_file_name.substr(0,str_file_name.find('.'));
    fprintf(of,"o %s\n",file_base_name.c_str());

    //====================================== Develop node mapping and write them
    std::vector<int> node_mapping(GetGlobalVertexCount(), -1);

    int node_counter=0;
    for (auto node : nodes_set)
    {
      node_counter++;
      int node_g_index = node;
      node_mapping[node_g_index] = node_counter;

      chi_mesh::Vertex cur_v = vertices[node_g_index];

      fprintf(of,"v %9.6f %9.6f %9.6f\n",cur_v.x,cur_v.y,cur_v.z);
    }

    //====================================== Write face normals
    for (const auto& face : faces_to_export)
    {
      fprintf(of,"vn %.4f %.4f %.4f\n",
              face.normal.x,
              face.normal.y,
              face.normal.z);
    }

    //====================================== Write faces
    int normal_counter=0;
    for (const auto& face : faces_to_export)
    {
      normal_counter++;
      fprintf(of,"f");

      for (auto v_g_index : face.vertex_ids)
        fprintf(of," %d//%d",node_mapping[v_g_index],normal_counter);

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
      str_file_name.substr(0,str_file_name.find('.'));

    if (chi_physics_handler.material_stack.empty())
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

      if (of == nullptr)
      {
        chi_log.Log(LOG_ALLERROR) << "Could not open file: "
                                  << mat_file_name;
        exit(EXIT_FAILURE);
      }

      //====================================== Develop list of faces and nodes
      std::set<int> nodes_set;
      std::vector<chi_mesh::CellFace> faces_to_export;
      for (const auto& cell : local_cells)
      {
        if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
        {
          if (cell.material_id != mat) continue;

          for (const auto& face : cell.faces)
          {
            int adjcell_glob_index = face.neighbor_id;

            if (adjcell_glob_index<0)
            {
              faces_to_export.push_back(face);

              for (auto vid : face.vertex_ids)
                nodes_set.insert(vid);
            }//if boundary
            else
            {
              auto& adj_cell = cells[adjcell_glob_index];

              if (adj_cell.material_id != mat)
              {
                faces_to_export.push_back(face);

                for (auto vid : face.vertex_ids)
                  nodes_set.insert(vid);
              }//if material missmatch
            }//if neigbor cell
          }//for face
        }//if polyhedron
      }//for local cell

      //====================================== Write header
      fprintf(of,"# Exported mesh file from Extrusion script\n");
      fprintf(of,"o %s\n",mat_base_name.c_str());

      //====================================== Develop node mapping and write them
      std::vector<int> node_mapping(GetGlobalVertexCount(), -1);

      int node_counter=0;
      for (auto node : nodes_set)
      {
        node_counter++;
        int node_g_index = node;
        node_mapping[node_g_index] = node_counter;

        chi_mesh::Vertex cur_v = vertices[node_g_index];

        fprintf(of,"v %9.6f %9.6f %9.6f\n",cur_v.x,cur_v.y,cur_v.z);
      }

      //====================================== Write face normals
      for (const auto& face : faces_to_export)
      {
        fprintf(of,"vn %.4f %.4f %.4f\n",
                face.normal.x,
                face.normal.y,
                face.normal.z);
      }

      //====================================== Write faces
      int normal_counter=0;
      for (const auto& face : faces_to_export)
      {
        normal_counter++;
        fprintf(of,"f");

        for (auto v_g_index : face.vertex_ids)
          fprintf(of," %d//%d",node_mapping[v_g_index],normal_counter);

        fprintf(of,"\n");
      }


      fclose(of);

      chi_log.Log(LOG_0)
        << "Exported Material Volume mesh to "
        << mat_file_name;
    }//for mat
  }//if per material

}