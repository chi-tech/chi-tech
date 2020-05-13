#include "surfacemesher.h"
#include "../MeshHandler/chi_meshhandler.h"
#include "../Region/chi_region.h"
#include "../Boundary/chi_boundary.h"
#include<iostream>

#include <chi_log.h>
extern ChiLog& chi_log;

//###################################################################
/**Virtual execute function. Meant to be overwritten.*/
void chi_mesh::SurfaceMesher::Execute()
{
  std::cout << "This is an empty mesher. Nothing to execute." << std::endl;
}


//###################################################################
/**Prints load balancing information based on a two D grid of XY-cuts.*/
void chi_mesh::SurfaceMesher::PrintLoadBalanceInfo()
{
  //================================================== Get the current handler
  chi_mesh::MeshHandler* mesh_handler = chi_mesh::GetCurrentHandler();

  //================================================== Loop over all regions
  std::vector<chi_mesh::Region*>::iterator region_iter;
  for (region_iter = mesh_handler->region_stack.begin();
       region_iter != mesh_handler->region_stack.end();
       region_iter++)
  {
    chi_mesh::Region* region = *region_iter;
    //====================================== Loop over boundaries
    std::vector<chi_mesh::Boundary*>::iterator bndry;
    for (bndry  = region->boundaries.begin();
         bndry != region->boundaries.end();
         bndry++)
    {
      if ((*bndry)->initial_mesh_continuum.surface_mesh != nullptr)
      {
        chi_mesh::SurfaceMesh* mesh =
          (*bndry)->initial_mesh_continuum.surface_mesh;

        int num_x_sets = xcuts.size()+1;
        int num_y_sets = ycuts.size()+1;
        int tot_num_sets = num_x_sets*num_y_sets;
        chi_log.Log(LOG_0VERBOSE_1) << "Total xy cell-sets: " << tot_num_sets;

        std::vector<int> set_counts(tot_num_sets,0);

        //================================== Triangle faces
        for (int f=0; f<mesh->faces.size(); f++)
        {
          //Compute centroid
          int vi0 = mesh->faces[f].v_index[0];
          int vi1 = mesh->faces[f].v_index[1];
          int vi2 = mesh->faces[f].v_index[2];

          chi_mesh::Vector3 v0 = mesh->vertices[vi0];
          chi_mesh::Vector3 v1 = mesh->vertices[vi1];
          chi_mesh::Vector3 v2 = mesh->vertices[vi2];

          chi_mesh::Vector3 c = (v0 + v1 + v2) / 3.0;

          //Compute partition
          int x=-1;
          for (int i=0; i<xcuts.size(); i++)
          {
            if (c.x <= xcuts[i])
            {x = i;break;}
          }
          if (x == -1) x=num_x_sets-1;

          int y=-1;
          for (int j=0; j<ycuts.size(); j++)
          {
            if (c.y <= ycuts[j])
            {y = j;break;}
          }
          if (y == -1) y=num_y_sets-1;


          int p = num_x_sets*y + x;

          set_counts[p] += 1;
        }//for tri faces

        //================================== Polygon faces
        for (int f=0; f<mesh->poly_faces.size(); f++)
        {
          //Compute centroid
          chi_mesh::Vector3 c;
          for (int v=0; v<mesh->poly_faces[f]->v_indices.size(); v++)
          {
            int vi = mesh->poly_faces[f]->v_indices[v];
            c = c + mesh->vertices[vi];
          }
          c = c/mesh->poly_faces[f]->v_indices.size();

          //Compute partition
          int x=-1;
          for (int i=0; i<xcuts.size(); i++)
          {
            if (c.x <= xcuts[i])
            {x = i;break;}
          }
          if (x == -1) x=num_x_sets-1;

          int y=-1;
          for (int j=0; j<ycuts.size(); j++)
          {
            if (c.y <= ycuts[j])
            {y = j;break;}
          }
          if (y == -1) y=num_y_sets-1;


          int p = num_x_sets*y + x;

          set_counts[p] += 1;
        }//for tri faces
        //======================================== Compute avg and max
        int total=0;
        int max=0;
        chi_log.Log(LOG_0)
        << "XY load balancing information:\n"
        << "  Set#   X#   Y#   #Cells";
        for (int y=0;y<num_y_sets; y++)
        {
          for (int x=0; x<num_x_sets; x++)
          {
            int p = num_y_sets*y + x;

            total += set_counts[p];
            if (set_counts[p]>max) max = set_counts[p];

            chi_log.Log(LOG_0)
              << std::right
              << std::setw(6) << p
              << std::setw(5) << x
              << std::setw(5) << y << " "
              << std::setw(6) << set_counts[p];
          }
        }
        double avg = total/set_counts.size();

        //======================================== Log info
        chi_log.Log(LOG_0)
        << "XY Load balance factor: "
        << std::setprecision(3) << std::fixed << max/avg;
      }//if surf mesh
    }//for bndry

  }//for region



}