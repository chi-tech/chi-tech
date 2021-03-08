#include "volmesher_extruder.h"
#include "ChiMesh/Cell/cell_polygon.h"
#include "ChiMesh/Cell/cell_polyhedron.h"

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Extrude template cells into polygons.*/
void chi_mesh::VolumeMesherExtruder::
  ExtrudeCells(chi_mesh::MeshContinuum& template_grid,
               chi_mesh::MeshContinuum& grid)
{
  //================================================== Start extrusion
  int num_global_cells = 0;
  for (int iz=0; iz<(vertex_layers.size()-1); iz++)
  {

    for (int tc=0; tc < template_grid.local_cells.size(); tc++)
    {
      //========================================= Get template cell
      if (template_grid.local_cells[tc].Type() !=
          chi_mesh::CellType::POLYGON)
      {
        chi_log.Log(LOG_ALLERROR) << "Extruder: Template cell error.";
        exit(EXIT_FAILURE);
      }
      auto template_cell = (chi_mesh::CellPolygon&)(template_grid.local_cells[tc]);


      //========================================= Precompute centroid
      auto centroid_precompd = ComputeTemplateCell3DCentroid(
        template_cell, template_grid, iz, iz + 1);

      //========================================= Create a template empty
      //                                          cell
      auto tcell = new chi_mesh::Cell(chi_mesh::CellType::GHOST);
      tcell->centroid = centroid_precompd;

      //========================================= Get the partition id
      tcell->partition_id =
        GetCellPartitionIDFromCentroid(tcell->centroid);

      //###################### NOT A LOCAL CELL ############################
      if ((tcell->partition_id != chi_mpi.location_id) and
          (!options.mesh_global))
      {
        tcell->global_id = num_global_cells;
        ++num_global_cells;

        bool is_neighbor_to_partition = IsTemplateCellNeighborToThisPartition(
          template_cell, template_grid, iz, tc);

        if (not is_neighbor_to_partition)
          delete tcell;
        else
        {
          grid.cells.push_back(tcell);
        }
      }
        //####################### LOCAL CELL ###############################
      else
      {
        //========================================= Create polyhedron
        auto cell = new chi_mesh::CellPolyhedron;
        cell->partition_id = tcell->partition_id;
        delete tcell;

        //========================================= Populate cell v-indices
        for (auto tc_vid : template_cell.vertex_ids)
          cell->vertex_ids.push_back(tc_vid + iz*node_z_index_incr);

        for (auto tc_vid : template_cell.vertex_ids)
          cell->vertex_ids.push_back(tc_vid + (iz+1)*node_z_index_incr);

        cell->centroid = centroid_precompd;
        cell->vertex_ids.shrink_to_fit();

        //========================================= Create side faces
        for (auto& face : template_cell.faces)
        {
          chi_mesh::CellFace newFace;

          newFace.vertex_ids.resize(4,-1);
          newFace.vertex_ids[0] = face.vertex_ids[0] + iz*node_z_index_incr;
          newFace.vertex_ids[1] = face.vertex_ids[1] + iz*node_z_index_incr;
          newFace.vertex_ids[2] = face.vertex_ids[1] + (iz+1)*node_z_index_incr;
          newFace.vertex_ids[3] = face.vertex_ids[0] + (iz+1)*node_z_index_incr;

          //Compute centroid
          chi_mesh::Vertex& v0 = *grid.vertices[newFace.vertex_ids[0]];
          chi_mesh::Vertex& v1 = *grid.vertices[newFace.vertex_ids[1]];
          chi_mesh::Vertex& v2 = *grid.vertices[newFace.vertex_ids[2]];
          chi_mesh::Vertex& v3 = *grid.vertices[newFace.vertex_ids[3]];

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
                                  iz*((int)template_grid.local_cells.size());
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
        newFace.vertex_ids.reserve(template_cell.vertex_ids.size());
        for (int tv=((int)(template_cell.vertex_ids.size())-1); tv>=0; tv--)
        {
          newFace.vertex_ids.push_back(template_cell.vertex_ids[tv]
                                       + iz*node_z_index_incr);
          chi_mesh::Vertex v = *grid.vertices[newFace.vertex_ids.back()];
          vfc = vfc + v;
        }

        //Compute centroid
        vfc = vfc/template_cell.vertex_ids.size();
        newFace.centroid = vfc;

        //Compute normal
        va = *grid.vertices[newFace.vertex_ids[0]] - vfc;
        vb = *grid.vertices[newFace.vertex_ids[1]] - vfc;

        vn = va.Cross(vb);
        newFace.normal = vn/vn.Norm();

        //Set neighbor
        if (iz==0)
        {
          newFace.neighbor_id = 0;
          newFace.has_neighbor = false;
        }
        else
        {
          newFace.neighbor_id = tc + (iz-1)*(int)(template_grid.local_cells.size());
          newFace.has_neighbor = true;
        }

        cell->faces.push_back(newFace);

        //=============================== Top face
        newFace = chi_mesh::CellFace();
        //Vertices
        vfc = chi_mesh::Vertex(0.0,0.0,0.0);
        newFace.vertex_ids.reserve(template_cell.vertex_ids.size());
        for (auto tc_vid : template_cell.vertex_ids)
        {
          newFace.vertex_ids.push_back(tc_vid + (iz+1)*node_z_index_incr);
          chi_mesh::Vertex v = *grid.vertices[newFace.vertex_ids.back()];
          vfc = vfc + v;
        }

        //Compute centroid
        vfc = vfc/template_cell.vertex_ids.size();
        newFace.centroid = vfc;

        //Compute normal
        va = *grid.vertices[newFace.vertex_ids[0]] - vfc;
        vb = *grid.vertices[newFace.vertex_ids[1]] - vfc;

        vn = va.Cross(vb);
        newFace.normal = vn/vn.Norm();

        //Set neighbor
        if (iz==(vertex_layers.size()-2))
        {
          newFace.neighbor_id = 0;
          newFace.has_neighbor = false;
        }
        else
        {
          newFace.neighbor_id = tc + (iz+1)*(int)(template_grid.local_cells.size());
          newFace.has_neighbor = true;
        }

        cell->faces.push_back(newFace);

        cell->global_id = num_global_cells;
        grid.cells.push_back(cell); ++num_global_cells;

      } //if cell is local
      //################### END OF LOCAL CELL ##############################

    }//for template cell

  }//for iz


}

