#include "volmesher_extruder.h"
#include "../../MeshHandler/chi_meshhandler.h"
#include "../../SurfaceMesher/surfacemesher.h"

#include <chi_log.h>
extern ChiLog chi_log;

#include <chi_mpi.h>
extern ChiMPI chi_mpi;

//###################################################################
/** Creates actual z-levels for the input layer specification.*/
void chi_mesh::VolumeMesherExtruder::SetupLayers(int default_layer_count)
{
  //================================================== Create default layers if no
  //                                                   input layers are provided
  if (input_layers.empty())
  {
    chi_log.Log(LOG_0WARNING)
      << "VolumeMesherExtruder: No extrusion layers have been specified. "
      << "A default single layer will be used with height 1.0 and a single "
      << "subdivision.";
    double dz = 1.0/default_layer_count;
    for (int k=0; k<=default_layer_count; k++)
    {
      vertex_layers.push_back(k*dz);
    }
  }
  else
  {
    double last_z=0.0;
    vertex_layers.push_back(last_z);

    for (int ell=0; ell<input_layers.size(); ell++)
    {
      double dz = input_layers[ell]->height/input_layers[ell]->sub_divisions;

      for (int k=0; k<input_layers[ell]->sub_divisions; k++)
      {
        last_z += dz;
        vertex_layers.push_back(last_z);
      }
    }
  }

  chi_log.Log(LOG_0)
    << "VolumeMesherExtruder: Total number of cell layers is "
    << vertex_layers.size()-1;
}

//###################################################################
/** Creates nodes that are owned locally.*/
void chi_mesh::VolumeMesherExtruder::
  CreateLocalAndBoundaryNodes(chi_mesh::MeshContinuum *template_continuum,
                              chi_mesh::MeshContinuum *vol_continuum)
{
  //================================================== Get current handler
  chi_mesh::MeshHandler* handler = chi_mesh::GetCurrentHandler();
  chi_mesh::SurfaceMesher* surf_mesher = handler->surface_mesher;

  //================================================== For each layer
  std::set<int> local_vert_ids;
  for (int iz=0; iz<(vertex_layers.size()-1); iz++)
  {
    for (int tc=0; tc<template_continuum->local_cells.size(); tc++)
    {
      //========================================= Get template cell
      if (template_continuum->local_cells[tc].Type() !=
          chi_mesh::CellType::POLYGON)
      {
        chi_log.Log(LOG_ALLERROR)
          << "Extruder::CreateLocalAndBoundaryNodes: Template cell error.";
        exit(EXIT_FAILURE);
      }
      auto template_cell = (chi_mesh::CellPolygon*)(&template_continuum->local_cells[tc]);

      //========================================= Precompute centroid
      auto centroid_precompd = ComputeTemplateCell3DCentroid(
        template_cell, template_continuum, iz, iz+1);

      //========================================= Get the partition id
      int tcell_partition_id =
        GetCellPartitionIDFromCentroid(centroid_precompd, surf_mesher);

      //###################### NOT A LOCAL CELL ############################
      if (tcell_partition_id != chi_mpi.location_id)
      {
        bool is_neighbor_to_partition = IsTemplateCellNeighborToThisPartition(
          template_cell, template_continuum, surf_mesher, iz, tc);

        if (is_neighbor_to_partition)
        {
          for (auto tc_vid : template_cell->vertex_ids)
            local_vert_ids.insert(tc_vid + iz*node_z_index_incr);

          for (auto tc_vid : template_cell->vertex_ids)
            local_vert_ids.insert(tc_vid + (iz+1)*node_z_index_incr);
        }
      }
      //####################### LOCAL CELL ###############################
      else
      {
        for (auto tc_vid : template_cell->vertex_ids)
          local_vert_ids.insert(tc_vid + iz*node_z_index_incr);

        for (auto tc_vid : template_cell->vertex_ids)
          local_vert_ids.insert(tc_vid + (iz+1)*node_z_index_incr);
      }
    }//for template cell
  }//for layer

  //============================================= Now add all nodes
  //                                              that are local or neighboring
  for (auto vert : vol_continuum->vertices) delete vert;
  vol_continuum->vertices.clear();
//  for (int iz=0; iz<vertex_layers.size(); iz++)
  for (auto layer_z_level : vertex_layers)
  {
    for (auto vertex : template_continuum->vertices)
    {
      int new_vert_index = vol_continuum->vertices.size();

      auto local_index = local_vert_ids.find(new_vert_index);

      if (local_index != local_vert_ids.end())
      {
        auto node = new chi_mesh::Node(*vertex);
        node->z = layer_z_level;

        vol_continuum->vertices.push_back(node);
      }
      else
        vol_continuum->vertices.push_back(nullptr);


    }
  }
}