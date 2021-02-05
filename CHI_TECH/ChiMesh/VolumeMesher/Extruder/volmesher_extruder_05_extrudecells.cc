#include "volmesher_extruder.h"
#include "ChiMesh/Cell/cell_polygon.h"
#include "ChiMesh/Cell/cell_polyhedron.h"
#include "../../MeshHandler/chi_meshhandler.h"
#include "../../SurfaceMesher/surfacemesher.h"

#include <chi_mpi.h>
extern ChiMPI& chi_mpi;

#include <chi_log.h>
extern ChiLog& chi_log;

//###################################################################
/**Computes the partition id of a template cell's projection onto 3D.*/
chi_mesh::Vector3
chi_mesh::VolumeMesherExtruder::ComputeTemplateCell3DCentroid(
  chi_mesh::CellPolygon *n_template_cell,
  chi_mesh::MeshContinuumPtr template_continuum,
  int z_level_begin,int z_level_end)
{
  chi_mesh::Vector3 n_centroid_precompd;
  size_t tc_num_verts = n_template_cell->vertex_ids.size();

   for (auto tc_vid : n_template_cell->vertex_ids)
  {
    auto temp_vert = *template_continuum->vertices[tc_vid];
    temp_vert.z = vertex_layers[z_level_begin];
    n_centroid_precompd = n_centroid_precompd + temp_vert;
  }

  for (auto tc_vid : n_template_cell->vertex_ids)
  {
    auto temp_vert = *template_continuum->vertices[tc_vid];
    temp_vert.z = vertex_layers[z_level_end];
    n_centroid_precompd = n_centroid_precompd + temp_vert;
  }
  n_centroid_precompd = n_centroid_precompd/(2*tc_num_verts);

  return n_centroid_precompd;
}

//###################################################################
/**Computes a cell's partition id based on a centroid.*/
int chi_mesh::VolumeMesherExtruder::
  GetCellPartitionIDFromCentroid(chi_mesh::Vector3& centroid,
                                 chi_mesh::SurfaceMesher* surf_mesher)
{
  int px = surf_mesher->partitioning_x;
  int py = surf_mesher->partitioning_y;

  chi_mesh::Cell n_gcell(chi_mesh::CellType::GHOST);
  n_gcell.centroid = centroid;

  auto xyz_partition_indices = GetCellXYZPartitionID(&n_gcell);

  int nxi = std::get<0>(xyz_partition_indices);
  int nyi = std::get<1>(xyz_partition_indices);
  int nzi = std::get<2>(xyz_partition_indices);

  return nzi*px*py + nyi*px + nxi;
}

//###################################################################
/**Determines if a template cell is neighbor to the current partition.*/
bool chi_mesh::VolumeMesherExtruder::
  IsTemplateCellNeighborToThisPartition(
    chi_mesh::CellPolygon *template_cell,
    chi_mesh::MeshContinuumPtr template_continuum,
    chi_mesh::SurfaceMesher *surf_mesher,
    int z_level,int tc_index)
{
  int iz = z_level;
  int tc = tc_index;

  //========================= Loop over template cell neighbors
  //                          for side neighbors
  bool is_neighbor_to_partition = false;
  for (auto& tc_face : template_cell->faces)
  {
    if (tc_face.has_neighbor)
    {
      auto n_template_cell = (chi_mesh::CellPolygon*)(
        &template_continuum->local_cells[tc_face.neighbor_id]);

      auto n_centroid_precompd = ComputeTemplateCell3DCentroid(
        n_template_cell, template_continuum, iz, iz+1);

      int n_gcell_partition_id =
        GetCellPartitionIDFromCentroid(n_centroid_precompd,surf_mesher);

      if (n_gcell_partition_id == chi_mpi.location_id)
      {
        is_neighbor_to_partition = true;
        break;
      }
    }//if neighbor not border
  }//for neighbors

  //========================= Now look at bottom neighbor
  //Bottom face
  if (iz != 0)
  {
    auto n_template_cell = (chi_mesh::CellPolygon*)
      (&template_continuum->local_cells[tc]);

    auto n_centroid_precompd = ComputeTemplateCell3DCentroid(
      n_template_cell, template_continuum, iz-1, iz);

    int n_gcell_partition_id =
      GetCellPartitionIDFromCentroid(n_centroid_precompd,surf_mesher);

    if (n_gcell_partition_id == chi_mpi.location_id)
      is_neighbor_to_partition = true;
  }//if neighbor not border

  //========================= Now look at top neighbor
  //Top Face
  if (iz != (vertex_layers.size()-2))
  {
    auto n_template_cell = (chi_mesh::CellPolygon*)
      (&template_continuum->local_cells[tc]);

    auto n_centroid_precompd = ComputeTemplateCell3DCentroid(
      n_template_cell, template_continuum, iz+1, iz+2);

    int n_gcell_partition_id =
      GetCellPartitionIDFromCentroid(n_centroid_precompd,surf_mesher);

    if (n_gcell_partition_id == chi_mpi.location_id)
      is_neighbor_to_partition = true;
  }//if neighbor not border

  return is_neighbor_to_partition;
}

//###################################################################
/**Extrude template cells into polygons.*/
void chi_mesh::VolumeMesherExtruder::
ExtrudeCells(chi_mesh::MeshContinuumPtr template_continuum,
             chi_mesh::MeshContinuumPtr vol_continuum)
{
  //================================================== Get current handler
  chi_mesh::MeshHandler* handler = chi_mesh::GetCurrentHandler();
  chi_mesh::SurfaceMesher* surf_mesher = handler->surface_mesher;

  //================================================== Start extrusion
  int num_global_cells = 0;
  for (int iz=0; iz<(vertex_layers.size()-1); iz++)
  {

    for (int tc=0; tc<template_continuum->local_cells.size(); tc++)
    {
      //========================================= Get template cell
      if (template_continuum->local_cells[tc].Type() !=
          chi_mesh::CellType::POLYGON)
      {
        chi_log.Log(LOG_ALLERROR) << "Extruder: Template cell error.";
        exit(EXIT_FAILURE);
      }
      auto template_cell = (chi_mesh::CellPolygon*)(&template_continuum->local_cells[tc]);


      //========================================= Precompute centroid
      auto centroid_precompd = ComputeTemplateCell3DCentroid(
        template_cell, template_continuum, iz, iz+1);

      //========================================= Create a template empty
      //                                          cell
      auto tcell = new chi_mesh::Cell(chi_mesh::CellType::GHOST);
      tcell->centroid = centroid_precompd;

      //========================================= Get the partition id
      tcell->partition_id =
        GetCellPartitionIDFromCentroid(tcell->centroid, surf_mesher);

      //###################### NOT A LOCAL CELL ############################
      if ((tcell->partition_id != chi_mpi.location_id) and
          (!options.mesh_global))
      {
        tcell->global_id = num_global_cells;
//        vol_continuum->cells.push_back(tcell); ++num_global_cells;
        ++num_global_cells;

        bool is_neighbor_to_partition = IsTemplateCellNeighborToThisPartition(
          template_cell, template_continuum, surf_mesher, iz, tc);

        if (not is_neighbor_to_partition)
        {
//          vol_continuum->cells[tcell->cell_global_id] = nullptr;
          delete tcell;
        }
        else
        {
          vol_continuum->cells.push_back(tcell);
        }
      }
        //####################### LOCAL CELL ###############################
      else
      {
        //========================================= Create polyhedron
        auto cell = new chi_mesh::CellPolyhedron;
//        cell->xyz_partition_indices = tcell->xyz_partition_indices;
        cell->partition_id = tcell->partition_id;
        delete tcell;

        //========================================= Populate cell v-indices
        for (auto tc_vid : template_cell->vertex_ids)
          cell->vertex_ids.push_back(tc_vid + iz*node_z_index_incr);

        for (auto tc_vid : template_cell->vertex_ids)
          cell->vertex_ids.push_back(tc_vid + (iz+1)*node_z_index_incr);

        cell->centroid = centroid_precompd;
        cell->vertex_ids.shrink_to_fit();

        //========================================= Create side faces
        for (auto& face : template_cell->faces)
        {
          chi_mesh::CellFace newFace;

          newFace.vertex_ids.resize(4,-1);
          newFace.vertex_ids[0] = face.vertex_ids[0] + iz*node_z_index_incr;
          newFace.vertex_ids[1] = face.vertex_ids[1] + iz*node_z_index_incr;
          newFace.vertex_ids[2] = face.vertex_ids[1] + (iz+1)*node_z_index_incr;
          newFace.vertex_ids[3] = face.vertex_ids[0] + (iz+1)*node_z_index_incr;

          //Compute centroid
          chi_mesh::Vertex& v0 = *vol_continuum->vertices[newFace.vertex_ids[0]];
          chi_mesh::Vertex& v1 = *vol_continuum->vertices[newFace.vertex_ids[1]];
          chi_mesh::Vertex& v2 = *vol_continuum->vertices[newFace.vertex_ids[2]];
          chi_mesh::Vertex& v3 = *vol_continuum->vertices[newFace.vertex_ids[3]];

          chi_mesh::Vertex vfc = (v0+v1+v2+v3)/4.0;
          newFace.centroid = vfc;

          //Compute normal
          chi_mesh::Vector3 va = v0 - vfc;
          chi_mesh::Vector3 vb = v1 - vfc;

          chi_mesh::Vector3 vn = va.Cross(vb);

          newFace.normal = (vn/vn.Norm());

          //Set neighbor
          //The side connections have the same connections as the
          //template cell + the iz specifiers of the layer.
          if (face.has_neighbor)
          {
            newFace.neighbor_id = face.neighbor_id +
                                  iz*((int)template_continuum->local_cells.size());
            newFace.has_neighbor = true;
          }
          else
            newFace.neighbor_id = face.neighbor_id;

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
          chi_mesh::Vertex v = *vol_continuum->vertices[newFace.vertex_ids.back()];
          vfc = vfc + v;
        }

        //Compute centroid
        vfc = vfc/template_cell->vertex_ids.size();
        newFace.centroid = vfc;

        //Compute normal
        va = *vol_continuum->vertices[newFace.vertex_ids[0]] - vfc;
        vb = *vol_continuum->vertices[newFace.vertex_ids[1]] - vfc;

        vn = va.Cross(vb);
        newFace.normal = vn/vn.Norm();

        //Set neighbor
        if (iz==0)
          newFace.neighbor_id = bot_boundary_index;
        else
        {
          newFace.neighbor_id = tc + (iz-1)*(int)(template_continuum->local_cells.size());
          newFace.has_neighbor = true;
        }

        cell->faces.push_back(newFace);

        //=============================== Top face
        newFace = chi_mesh::CellFace();
        //Vertices
        vfc = chi_mesh::Vertex(0.0,0.0,0.0);
        newFace.vertex_ids.reserve(template_cell->vertex_ids.size());
        for (auto tc_vid : template_cell->vertex_ids)
        {
          newFace.vertex_ids.push_back(tc_vid + (iz+1)*node_z_index_incr);
          chi_mesh::Vertex v = *vol_continuum->vertices[newFace.vertex_ids.back()];
          vfc = vfc + v;
        }

        //Compute centroid
        vfc = vfc/template_cell->vertex_ids.size();
        newFace.centroid = vfc;

        //Compute normal
        va = *vol_continuum->vertices[newFace.vertex_ids[0]] - vfc;
        vb = *vol_continuum->vertices[newFace.vertex_ids[1]] - vfc;

        vn = va.Cross(vb);
        newFace.normal = vn/vn.Norm();

        //Set neighbor
        if (iz==(vertex_layers.size()-2))
          newFace.neighbor_id = top_boundary_index;
        else
        {
          newFace.neighbor_id = tc + (iz+1)*(int)(template_continuum->local_cells.size());
          newFace.has_neighbor = true;
        }

        cell->faces.push_back(newFace);

        cell->global_id = num_global_cells;
        vol_continuum->cells.push_back(cell); ++num_global_cells;

      } //if cell is local
      //################### END OF LOCAL CELL ##############################

    }//for template cell

  }//for iz


}

