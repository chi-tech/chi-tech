#include "volmesher_extruder.h"
#include "../../CHI_MESHCONTINUUM/chi_meshcontinuum.h"
#include "../../CHI_CELL/cell_polyhedron.h"
#include "../../CHI_CELL/cell_polygon.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../CHI_SURFACEMESHER/surfacemesher.h"
#include "../../../CHI_MPI/chi_mpi.h"

extern CHI_MPI chi_mpi;

//###################################################################
/**Extrude template cells into polygons.*/
void chi_mesh::VolumeMesherExtruder::
ExtrudeCells(chi_mesh::MeshContinuum *template_continuum,
             chi_mesh::MeshContinuum *vol_continuum)
{
  //================================================== Get current handler
  chi_mesh::MeshHandler* handler = chi_mesh::GetCurrentHandler();
  chi_mesh::SurfaceMesher* surf_mesher = handler->surface_mesher;

  //================================================== Start extrusion
  for (int iz=0; iz<(vertex_layers.size()-1); iz++)
  {

    for (int tc=0; tc<template_continuum->cells.size(); tc++)
    {
      //========================================= Get template cell
      chi_mesh::CellPolygon* template_cell =
        (chi_mesh::CellPolygon*)template_continuum->cells[tc];

      //========================================= Precompute centroid
      chi_mesh::Vector centroid_precompd;
      for (int tv=0; tv<template_cell->v_indices.size(); tv++)
      {
        int v_index = (template_cell->v_indices[tv]
                       + iz*node_z_index_incr);
        centroid_precompd = centroid_precompd + *vol_continuum->nodes[v_index];
      }
      for (int tv=0; tv<template_cell->v_indices.size(); tv++)
      {
        int v_index = (template_cell->v_indices[tv]
                       + (iz+1)*node_z_index_incr);
        centroid_precompd = centroid_precompd + *vol_continuum->nodes[(v_index)];
      }
      centroid_precompd = centroid_precompd/(2*template_cell->v_indices.size());

      //========================================= Create a template empty
      //                                          cell
      chi_mesh::Cell *tcell = new chi_mesh::Cell;
      tcell->centroid = centroid_precompd;

      //========================================= Get the partition id
      tcell->xyz_partition_indices = GetCellXYZPartitionID(tcell);
      int xi,yi,zi;
      xi = std::get<0>(tcell->xyz_partition_indices);
      yi = std::get<1>(tcell->xyz_partition_indices);
      zi = std::get<2>(tcell->xyz_partition_indices);
      int px,py;
      px = surf_mesher->partitioning_x;
      py = surf_mesher->partitioning_y;
      tcell->partition_id = zi*px*py + yi*px + xi;
      //printf("xi,yi,zi=%d,%d,%d\n", xi,yi,zi);

      //###################### NOT A LOCAL CELL ############################
      if (tcell->partition_id != chi_mpi.location_id)
      {
        tcell->cell_global_id = vol_continuum->cells.size();
        vol_continuum->cells.push_back(tcell);

      }
        //####################### LOCAL CELL ###############################
      else
      {
        delete [] tcell;
        //========================================= Create polyhedron
        chi_mesh::CellPolyhedron* cell = new chi_mesh::CellPolyhedron;
        cell->xyz_partition_indices = tcell->xyz_partition_indices;
        cell->partition_id = tcell->partition_id;

        //========================================= Compute centroid
        //                                          and set vindices
        for (int tv=0; tv<template_cell->v_indices.size(); tv++)
        {
          int v_index = (template_cell->v_indices[tv]
                         + iz*node_z_index_incr);
          cell->v_indices.push_back(v_index);
        }
        for (int tv=0; tv<template_cell->v_indices.size(); tv++)
        {
          int v_index = (template_cell->v_indices[tv]
                         + (iz+1)*node_z_index_incr);
          cell->v_indices.push_back(v_index);
        }
        cell->centroid = centroid_precompd;
        cell->v_indices.shrink_to_fit();

        //========================================= Create side faces
        for (int e=0; e<template_cell->edges.size(); e++)
        {
          int* template_edge_indices = template_cell->edges[e];

          chi_mesh::PolyFace* newFace = new chi_mesh::PolyFace;

          newFace->v_indices.push_back(template_edge_indices[0]
                                       +iz*node_z_index_incr);
          newFace->v_indices.push_back(template_edge_indices[1]
                                       +iz*node_z_index_incr);
          newFace->v_indices.push_back(template_edge_indices[1]
                                       +(iz+1)*node_z_index_incr);
          newFace->v_indices.push_back(template_edge_indices[0]
                                       +(iz+1)*node_z_index_incr);
          newFace->v_indices.shrink_to_fit();

          int* newEdge;

          newEdge = new int[4];
          newEdge[0] = newFace->v_indices[0];
          newEdge[1] = newFace->v_indices[1];
          newEdge[2] = -1;
          newEdge[3] = -1;
          newFace->edges.push_back(newEdge);

          newEdge = new int[4];
          newEdge[0] = newFace->v_indices[1];
          newEdge[1] = newFace->v_indices[2];
          newEdge[2] = -1;
          newEdge[3] = -1;
          newFace->edges.push_back(newEdge);

          newEdge = new int[4];
          newEdge[0] = newFace->v_indices[2];
          newEdge[1] = newFace->v_indices[3];
          newEdge[2] = -1;
          newEdge[3] = -1;
          newFace->edges.push_back(newEdge);

          newEdge = new int[4];
          newEdge[0] = newFace->v_indices[3];
          newEdge[1] = newFace->v_indices[0];
          newEdge[2] = -1;
          newEdge[3] = -1;
          newFace->edges.push_back(newEdge);

          newFace->edges.shrink_to_fit();

          //Compute centroid and normal
          chi_mesh::Vertex v0 = *vol_continuum->nodes[newFace->v_indices[0]];
          chi_mesh::Vertex v1 = *vol_continuum->nodes[newFace->v_indices[1]];
          chi_mesh::Vertex v2 = *vol_continuum->nodes[newFace->v_indices[2]];
          chi_mesh::Vertex v3 = *vol_continuum->nodes[newFace->v_indices[3]];

          chi_mesh::Vertex vfc = (v0+v1+v2+v3)/4.0;
          newFace->face_centroid = vfc;

          chi_mesh::Vector va = v0 - vfc;
          chi_mesh::Vector vb = v1 - vfc;

          chi_mesh::Vector vn = va.Cross(vb);

          newFace->geometric_normal = (vn/vn.Norm());


//          if (chi_mpi.location_id == 0)
//          {
//            std::cout << "v0=" << v0.PrintS() << "\n";
//            std::cout << "v1=" << v1.PrintS() << "\n";
//            std::cout << "v2=" << v2.PrintS() << "\n";
//            std::cout << "v3=" << v3.PrintS() << "\n";
//            std::cout << "n" << cell->faces.size() << "=" << newFace->geometric_normal.PrintS() << "\n";
//          }


          //The side connections have the same connections as the
          //template cell + the iz specifiers of the layer.
          if (template_edge_indices[2]<0)
          {
            newFace->face_indices[0] = template_edge_indices[2];
            newFace->face_indices[1] = template_edge_indices[3];
          }
          else
          {
            newFace->face_indices[0] = template_edge_indices[2]
                                       + iz*((int)template_continuum->cells.size());
            newFace->face_indices[1] = template_edge_indices[3];
          }


          cell->faces.push_back(newFace);
        } //for side faces

        //========================================= Create top and bottom faces
        chi_mesh::PolyFace* newFace;

        chi_mesh::Vertex vfc;
        chi_mesh::Vertex va;
        chi_mesh::Vertex vb;
        chi_mesh::Vertex vn;

        //=============================== Bottom] face
        newFace = new chi_mesh::PolyFace;
        //Vertices
        vfc = chi_mesh::Vertex(0.0,0.0,0.0);
        for (int tv=((int)(template_cell->v_indices.size())-1); tv>=0; tv--)
        {
          newFace->v_indices.push_back(template_cell->v_indices[tv]
                                       + iz*node_z_index_incr);
          chi_mesh::Vertex v = *vol_continuum->nodes[newFace->v_indices.back()];
          vfc = vfc + v;

          if (iz==0)
          {
            newFace->face_indices[0] = -1*(bot_boundary_index+1);
            newFace->face_indices[1] = -1;
          }
          else
          {
            newFace->face_indices[0] = tc +
                               (iz-1)*(int)(template_continuum->cells.size());
            newFace->face_indices[1] = (int)(template_cell->edges.size())+2-1;
          }
        }
        vfc = vfc/template_cell->v_indices.size();
        newFace->face_centroid = vfc;
        va = *vol_continuum->nodes[newFace->v_indices[0]] - vfc;
        vb = *vol_continuum->nodes[newFace->v_indices[1]] - vfc;

        vn = va.Cross(vb);
        newFace->geometric_normal = vn/vn.Norm();

//        std::cout << "bn=" << newFace->geometric_normal.PrintS() << "\n";

        newFace->v_indices.shrink_to_fit();
        //Edges
        for (int v=0; v<newFace->v_indices.size(); v++)
        {
          int* newEdge;
          newEdge = new int[4];
          newEdge[0] = newFace->v_indices[v];
          if (v == (newFace->v_indices.size()-1))
          {
            newEdge[1] = newFace->v_indices[0];
          }
          else
          {
            newEdge[1] = newFace->v_indices[v+1];
          }
          newEdge[2] = -1;
          newEdge[3] = -1;
          newFace->edges.push_back(newEdge);
        }
        newFace->edges.shrink_to_fit();
        cell->faces.push_back(newFace);

        //=============================== Top face
        newFace = new chi_mesh::PolyFace;
        //Vertices
        vfc = chi_mesh::Vertex(0.0,0.0,0.0);
        for (int tv=0; tv<template_cell->v_indices.size(); tv++)
        {
          newFace->v_indices.push_back(template_cell->v_indices[tv]
                                       + (iz+1)*node_z_index_incr);
          chi_mesh::Vertex v = *vol_continuum->nodes[newFace->v_indices.back()];
          vfc = vfc + v;

          if (iz==(vertex_layers.size()-2))
          {
            newFace->face_indices[0] = -1*(top_boundary_index+1);
            newFace->face_indices[1] = -1;
          }
          else
          {
            newFace->face_indices[0] = tc +
                           (iz+1)*(int)(template_continuum->cells.size());
            newFace->face_indices[1] = (int)(template_cell->edges.size())+1-1;
          }
        }
        vfc = vfc/template_cell->v_indices.size();
        newFace->face_centroid = vfc;
        va = *vol_continuum->nodes[newFace->v_indices[0]] - vfc;
        vb = *vol_continuum->nodes[newFace->v_indices[1]] - vfc;

        vn = va.Cross(vb);
        newFace->geometric_normal = vn/vn.Norm();

//        std::cout << "tn=" << newFace->geometric_normal.PrintS() << "\n";

        newFace->v_indices.shrink_to_fit();
        //Edges
        for (int v=0; v<newFace->v_indices.size(); v++)
        {
          int* newEdge;
          newEdge = new int[4];
          newEdge[0] = newFace->v_indices[v];
          if (v == (newFace->v_indices.size()-1))
          {
            newEdge[1] = newFace->v_indices[0];
          }
          else
          {
            newEdge[1] = newFace->v_indices[v+1];
          }
          newEdge[2] = -1;
          newEdge[3] = -1;
          newFace->edges.push_back(newEdge);
        }
        newFace->edges.shrink_to_fit();
        cell->faces.push_back(newFace);


        cell->cell_global_id = vol_continuum->cells.size();
        vol_continuum->cells.push_back(cell);


      } //if cell is local
      //################### END OF LOCAL CELL ##############################

    }//for template cell

  }//for iz


}