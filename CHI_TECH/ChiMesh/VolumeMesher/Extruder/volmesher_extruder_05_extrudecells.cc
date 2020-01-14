#include "volmesher_extruder.h"
#include "../../MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell_polygon.h"
#include "ChiMesh/Cell/cell_polyhedron.h"
#include "../../MeshHandler/chi_meshhandler.h"
#include "../../SurfaceMesher/surfacemesher.h"

#include <chi_mpi.h>
extern ChiMPI chi_mpi;

#include <chi_log.h>
extern ChiLog chi_log;

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
      auto template_cell =
        dynamic_cast<chi_mesh::CellPolygon*>(template_continuum->cells[tc]);
//      chi_mesh::CellPolygon* template_cell =
//        (chi_mesh::CellPolygon*)template_continuum->cells[tc];
      if (template_cell == nullptr)
      {
        chi_log.Log(LOG_ALLERROR) << "Extruder: Template cell error.";
        exit(EXIT_FAILURE);
      }

      //========================================= Precompute centroid
      chi_mesh::Vector centroid_precompd;
      for (int tv=0; tv<template_cell->vertex_ids.size(); tv++)
      {
        int v_index = (template_cell->vertex_ids[tv]
                       + iz*node_z_index_incr);
        centroid_precompd = centroid_precompd + *vol_continuum->nodes[v_index];
      }
      for (int tv=0; tv<template_cell->vertex_ids.size(); tv++)
      {
        int v_index = (template_cell->vertex_ids[tv]
                       + (iz+1)*node_z_index_incr);
        centroid_precompd = centroid_precompd + *vol_continuum->nodes[(v_index)];
      }
      centroid_precompd = centroid_precompd/(2*template_cell->vertex_ids.size());

      //========================================= Create a template empty
      //                                          cell
      chi_mesh::Cell *tcell = new chi_mesh::Cell(chi_mesh::CellType::GHOST);
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
      if ((tcell->partition_id != chi_mpi.location_id) and
          (!options.mesh_global))
      {
        tcell->cell_global_id = vol_continuum->cells.size();
        vol_continuum->cells.push_back(tcell);

      }
        //####################### LOCAL CELL ###############################
      else
      {
        delete tcell;
        //========================================= Create polyhedron
        chi_mesh::CellPolyhedron* cell = new chi_mesh::CellPolyhedron;
        cell->xyz_partition_indices = tcell->xyz_partition_indices;
        cell->partition_id = tcell->partition_id;

        //========================================= Compute centroid
        //                                          and set vindices
        for (int tv=0; tv<template_cell->vertex_ids.size(); tv++)
        {
          int v_index = (template_cell->vertex_ids[tv]
                         + iz*node_z_index_incr);
          cell->vertex_ids.push_back(v_index);
        }
        for (int tv=0; tv<template_cell->vertex_ids.size(); tv++)
        {
          int v_index = (template_cell->vertex_ids[tv]
                         + (iz+1)*node_z_index_incr);
          cell->vertex_ids.push_back(v_index);
        }
        cell->centroid = centroid_precompd;
        cell->vertex_ids.shrink_to_fit();

        //========================================= Create side faces
        for (int e=0; e<template_cell->faces.size(); e++)
        {
          chi_mesh::CellFace& face = template_cell->faces[e];

          chi_mesh::CellFace newFace;

          newFace.vertex_ids.resize(4,-1);
          newFace.vertex_ids[0] = face.vertex_ids[0] + iz*node_z_index_incr;
          newFace.vertex_ids[1] = face.vertex_ids[1] + iz*node_z_index_incr;
          newFace.vertex_ids[2] = face.vertex_ids[1] + (iz+1)*node_z_index_incr;
          newFace.vertex_ids[3] = face.vertex_ids[0] + (iz+1)*node_z_index_incr;

          //Compute centroid
          chi_mesh::Vertex v0 = *vol_continuum->nodes[newFace.vertex_ids[0]];
          chi_mesh::Vertex v1 = *vol_continuum->nodes[newFace.vertex_ids[1]];
          chi_mesh::Vertex v2 = *vol_continuum->nodes[newFace.vertex_ids[2]];
          chi_mesh::Vertex v3 = *vol_continuum->nodes[newFace.vertex_ids[3]];

          chi_mesh::Vertex vfc = (v0+v1+v2+v3)/4.0;
          newFace.centroid = vfc;

          //Compute normal
          chi_mesh::Vector va = v0 - vfc;
          chi_mesh::Vector vb = v1 - vfc;

          chi_mesh::Vector vn = va.Cross(vb);

          newFace.normal = (vn/vn.Norm());

          //Set neighbor
          //The side connections have the same connections as the
          //template cell + the iz specifiers of the layer.
          if (face.neighbor >= 0)
            newFace.neighbor = face.neighbor +
                               iz*((int)template_continuum->cells.size());
          else
            newFace.neighbor = face.neighbor;

          cell->faces.push_back(newFace);
        } //for side faces

        //========================================= Create top and bottom faces
        chi_mesh::CellFace newFace;

        chi_mesh::Vertex vfc;
        chi_mesh::Vertex va;
        chi_mesh::Vertex vb;
        chi_mesh::Vertex vn;

        //=============================== Bottom face
        newFace = chi_mesh::CellFace();
        //Vertices
        vfc = chi_mesh::Vertex(0.0,0.0,0.0);
        newFace.vertex_ids.reserve(template_cell->vertex_ids.size());
        for (int tv=((int)(template_cell->vertex_ids.size())-1); tv>=0; tv--)
        {
          newFace.vertex_ids.push_back(template_cell->vertex_ids[tv]
                                       + iz*node_z_index_incr);
          chi_mesh::Vertex v = *vol_continuum->nodes[newFace.vertex_ids.back()];
          vfc = vfc + v;
        }

        //Compute centroid
        vfc = vfc/template_cell->vertex_ids.size();
        newFace.centroid = vfc;

        //Compute normal
        va = *vol_continuum->nodes[newFace.vertex_ids[0]] - vfc;
        vb = *vol_continuum->nodes[newFace.vertex_ids[1]] - vfc;

        vn = va.Cross(vb);
        newFace.normal = vn/vn.Norm();

        //Set neighbor
        if (iz==0)
          newFace.neighbor = -1*(bot_boundary_index+1);
        else
          newFace.neighbor = tc + (iz-1)*(int)(template_continuum->cells.size());

        cell->faces.push_back(newFace);

        //=============================== Top face
        newFace = chi_mesh::CellFace();
        //Vertices
        vfc = chi_mesh::Vertex(0.0,0.0,0.0);
        newFace.vertex_ids.reserve(template_cell->vertex_ids.size());
        for (int tv=0; tv<template_cell->vertex_ids.size(); tv++)
        {
          newFace.vertex_ids.push_back(template_cell->vertex_ids[tv]
                                       + (iz+1)*node_z_index_incr);
          chi_mesh::Vertex v = *vol_continuum->nodes[newFace.vertex_ids.back()];
          vfc = vfc + v;
        }

        //Compute centroid
        vfc = vfc/template_cell->vertex_ids.size();
        newFace.centroid = vfc;

        //Compute normal
        va = *vol_continuum->nodes[newFace.vertex_ids[0]] - vfc;
        vb = *vol_continuum->nodes[newFace.vertex_ids[1]] - vfc;

        vn = va.Cross(vb);
        newFace.normal = vn/vn.Norm();

        //Set neighbor
        if (iz==(vertex_layers.size()-2))
          newFace.neighbor = -1*(top_boundary_index+1);
        else
          newFace.neighbor = tc + (iz+1)*(int)(template_continuum->cells.size());

        cell->faces.push_back(newFace);

        cell->cell_global_id = vol_continuum->cells.size();
        vol_continuum->cells.push_back(cell);

      } //if cell is local
      //################### END OF LOCAL CELL ##############################

    }//for template cell

  }//for iz


}

