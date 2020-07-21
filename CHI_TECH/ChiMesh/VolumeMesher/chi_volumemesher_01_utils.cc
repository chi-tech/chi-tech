#include "chi_volumemesher.h"
#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include <ChiMesh/Region/chi_region.h>
#include <ChiMesh/SurfaceMesh/chi_surfacemesh.h>
#include <ChiMesh/SurfaceMesher/Predefined/surfmesher_predefined.h>
#include <ChiMesh/VolumeMesher/Predefined2D/volmesher_predefined2d.h>
#include <ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h>
#include <ChiMesh/VolumeMesher/Predefined3D/volmesher_predefined3d.h>
#include "Linemesh1D/volmesher_linemesh1d.h"

#include <ChiMesh/Cell/cell_polygon.h>
#include "../MeshHandler/chi_meshhandler.h"
#include "../../ChiMPI/chi_mpi.h"
#include "../LogicalVolume/chi_mesh_logicalvolume.h"

#include <chi_log.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

#include <ChiTimer/chi_timer.h>
extern ChiTimer chi_program_timer;

//###################################################################
/**Creates 2D polygon cells for each face of a surface mesh.*/
void chi_mesh::VolumeMesher::
  AddContinuumToRegion(chi_mesh::MeshContinuum *grid, chi_mesh::Region& region)
{
  region.volume_mesh_continua.push_back(grid);
}

//###################################################################
/**Creates 2D polygon cells for each face of a surface mesh.*/
void chi_mesh::VolumeMesher::
CreatePolygonCells(chi_mesh::SurfaceMesh *surface_mesh,
                   chi_mesh::MeshContinuum *vol_continuum,
                   bool delete_surface_mesh_elements,
                   bool force_local)
{
  //============================================= Get current mesh handler
  chi_mesh::MeshHandler* handler = chi_mesh::GetCurrentHandler();

  //============================================= Copy nodes
  for (auto& vertex : surface_mesh->vertices)
    vol_continuum->vertices.push_back(new chi_mesh::Node(vertex));

  //============================================= Delete nodes
  if (delete_surface_mesh_elements)
    surface_mesh->vertices = std::move(std::vector<chi_mesh::Vertex>(0));

  //============================================= Process faces
  unsigned int num_cells = 0;
  for (auto& face : surface_mesh->faces)
  {
    auto cell = new chi_mesh::CellPolygon;

    for (int k=0;k<3;k++)
    {
      cell->vertex_ids.push_back(face.v_index[k]);

      chi_mesh::CellFace new_face;

      new_face.vertex_ids.push_back(face.e_index[k][0]);
      new_face.vertex_ids.push_back(face.e_index[k][1]);


      chi_mesh::Vertex& v0 = *vol_continuum->vertices[face.e_index[k][0]];
      chi_mesh::Vertex& v1 = *vol_continuum->vertices[face.e_index[k][1]];
      new_face.centroid = v0*0.5 + v1*0.5;

      chi_mesh::Vector3 vk = chi_mesh::Vector3(0.0, 0.0, 1.0);

      chi_mesh::Vector3 va = v1 - v0;
      chi_mesh::Vector3 vn = va.Cross(vk);
      vn = vn/vn.Norm();
      new_face.normal = vn;

      new_face.neighbor = face.e_index[k][2];

      cell->faces.push_back(new_face);

      cell->centroid = cell->centroid + surface_mesh->vertices[face.v_index[k]];
    }
    cell->centroid = cell->centroid/3;

    //====================================== Compute xy partition id
    cell->xy_partition_indices = GetCellXYPartitionID(cell);
    cell->partition_id = cell->xy_partition_indices.second*
                         handler->surface_mesher->partitioning_x +
                         cell->xy_partition_indices.first;

    if (force_local)
      cell->partition_id = chi_mpi.location_id;

    cell->global_id = num_cells;

    vol_continuum->cells.push_back(cell); ++num_cells;
  }

  for (auto face : surface_mesh->poly_faces)
  {
    auto cell = new chi_mesh::CellPolygon;

    //====================================== Copy vertices
    for (auto vid : face->v_indices)
    {
      cell->vertex_ids.push_back(vid);
      cell->centroid = cell->centroid +
                       *vol_continuum->vertices[vid];
    }
    cell->centroid = cell->centroid/cell->vertex_ids.size();

    //====================================== Copy edges
    for (auto src_side : face->edges)
    {
      chi_mesh::CellFace new_face;

      new_face.vertex_ids.push_back(src_side[0]);
      new_face.vertex_ids.push_back(src_side[1]);

      chi_mesh::Vertex& v0 = *vol_continuum->vertices[src_side[0]];
      chi_mesh::Vertex& v1 = *vol_continuum->vertices[src_side[1]];
      new_face.centroid = v0*0.5 + v1*0.5;
      chi_mesh::Vector3 vk = chi_mesh::Vector3(0.0, 0.0, 1.0);

      chi_mesh::Vector3 va = v1 - v0;
      chi_mesh::Vector3 vn = va.Cross(vk);
      vn = vn/vn.Norm();
      new_face.normal = vn;

      new_face.neighbor = src_side[2];

      cell->faces.push_back(new_face);
    }

    //====================================== Compute partition id
    cell->xy_partition_indices = GetCellXYPartitionID(cell);
    cell->partition_id = cell->xy_partition_indices.second*
                         handler->surface_mesher->partitioning_x +
                         cell->xy_partition_indices.first;

    if (force_local)
      cell->partition_id = chi_mpi.location_id;

    cell->global_id = num_cells;

    vol_continuum->cells.push_back(cell); ++num_cells;

    if (delete_surface_mesh_elements)
      delete face;
  }

  if (delete_surface_mesh_elements)
    surface_mesh->poly_faces.clear();

}

//###################################################################
/**Filters non-esential ghosts from the grid.*/
void chi_mesh::VolumeMesher::
  GridFilterGhosts(chi_mesh::MeshContinuum *in_grid,
                   chi_mesh::MeshContinuum *out_grid)
{
  chi_log.Log(LOG_0VERBOSE_1) << "Filtering ghosts.";
  //======================================== Copy vertices
  for (auto vertex : in_grid->vertices)
    out_grid->vertices.push_back(vertex);

  //======================================== Copy local cells
  for (auto& cell : in_grid->local_cells)
    out_grid->cells.push_back(&cell);



  //======================================== Copy ghost cells only
  //                                         if neighbor to current partition
  auto in_ghost_ids = in_grid->cells.GetGhostGlobalIDs();
  std::vector<Cell*> cells_to_delete;
  for (int ghost_id : in_ghost_ids)
  {
    auto ref_ghost_cell = in_grid->cells[ghost_id];

    bool is_neighbor_to_this_loc = false;
    for (auto& face : ref_ghost_cell->faces)
    {
      if (face.neighbor<0) continue;

      auto adj_cell = in_grid->cells[face.neighbor];
      if (adj_cell->partition_id == chi_mpi.location_id)
      {
        is_neighbor_to_this_loc = true;
        break;
      }
    }//for face

    if (is_neighbor_to_this_loc)
      out_grid->cells.push_back(ref_ghost_cell);
    else
      cells_to_delete.push_back(ref_ghost_cell);
  }//for cell


  //======================================== Deleting non-essential ghosts
  for (auto cell : cells_to_delete)
    delete cell;


  chi_log.Log(LOG_0VERBOSE_1) << "Done filtering ghosts.";
}




//###################################################################
/**Obtains the xy partition IDs of a cell.
 * Cell xy_partition ids are obtained from
 * the surface mesher.*/
std::pair<int,int> chi_mesh::VolumeMesher::
 GetCellXYPartitionID(chi_mesh::Cell *cell)
{
  std::pair<int,int> ij_id(0,0);
  bool found_partition = false;

  if (chi_mpi.process_count == 1){return ij_id;}

  //================================================== Get the current handler
  chi_mesh::MeshHandler*  mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::SurfaceMesher* surf_mesher = mesh_handler->surface_mesher;

  if  (typeid(*surf_mesher) == typeid(chi_mesh::SurfaceMesherPredefined))
  {
    //====================================== Sanity check on partitioning
    int num_x_subsets = surf_mesher->xcuts.size()+1;
    int num_y_subsets = surf_mesher->ycuts.size()+1;

    int x_remainder = num_x_subsets%surf_mesher->partitioning_x;
    int y_remainder = num_y_subsets%surf_mesher->partitioning_y;

    if (x_remainder != 0)
    {
      chi_log.Log(LOG_ALLERROR)
        << "When specifying x-partitioning, the number of grp_subsets in x "
           "needs to be divisible by the number of partitions in x.";
      exit(EXIT_FAILURE);
    }

    if (y_remainder != 0)
    {
      chi_log.Log(LOG_ALLERROR)
        << "When specifying y-partitioning, the number of grp_subsets in y "
           "needs to be divisible by the number of partitions in y.";
      exit(EXIT_FAILURE);
    }

    int subsets_per_partitionx = num_x_subsets/surf_mesher->partitioning_x;
    int subsets_per_partitiony = num_y_subsets/surf_mesher->partitioning_y;

//    chi_log.Log(LOG_0ERROR) << num_x_subsets;
//    chi_log.Log(LOG_0ERROR) << num_y_subsets;
//    chi_log.Log(LOG_0ERROR) << subsets_per_partitionx;
//    chi_log.Log(LOG_0ERROR) << subsets_per_partitiony;


    //====================================== Determine x-partition
    int x=-1;
    int xcount=-1;
    for (int i =  subsets_per_partitionx-1;
             i <  surf_mesher->xcuts.size();
             i += subsets_per_partitionx)
    {
      xcount++;
      if (cell->centroid.x <= surf_mesher->xcuts[i])
      {
        x = xcount;
        break;
      }
    }
    if (x<0)
    {
      x = surf_mesher->partitioning_x-1;
    }

    //====================================== Determine y-partition
    int y=-1;
    int ycount=-1;
    for (int i =  subsets_per_partitiony-1;
         i <  surf_mesher->ycuts.size();
         i += subsets_per_partitiony)
    {
      ycount++;
      if (cell->centroid.y <= surf_mesher->ycuts[i])
      {
        y = ycount;
        break;
      }
    }
    if (y<0)
    {
      y = surf_mesher->partitioning_y - 1;
    }

    //====================================== Set partitioning
    ij_id.first = x;
    ij_id.second= y;

  }//if typeid

  return ij_id;
}

//###################################################################
/**Obtains the xyz partition IDs of a cell.
 * Cell xy_partition ids are obtained from
 * the surface mesher. z id is obtained from the volume mesher.*/
std::tuple<int,int,int> chi_mesh::VolumeMesher::
  GetCellXYZPartitionID(chi_mesh::Cell *cell)
{
  std::tuple<int,int,int> ijk_id(0,0,0);
  bool found_partition = false;

  if (chi_mpi.process_count == 1){return ijk_id;}

  //================================================== Get ij indices
  std::pair<int,int> ij_id = GetCellXYPartitionID(cell);


  //================================================== Get the current handler
  chi_mesh::MeshHandler*  mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher* vol_mesher = mesh_handler->volume_mesher;

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& SLAB
  if (typeid(*vol_mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
  {
    auto line_mesher = (chi_mesh::VolumeMesherLinemesh1D*)vol_mesher;

    if (chi_mpi.process_count != options.partition_z and !options.mesh_global)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Number of process requested, " << options.partition_z
        << ", in PARTITION_Z does not match the amount of processes "
        << "available " << chi_mpi.process_count;
      exit(EXIT_FAILURE);
    }

    int cells_per_loc =
      ceil(line_mesher->num_slab_cells/(double)options.partition_z);

    int cur_loc = 0;
    for (int k=0; k<chi_mpi.process_count; k++)
    {
      if (cell->global_id < ((k + 1) * cells_per_loc))
      {
        cur_loc = k;
        found_partition = true;
        break;
      }
    }//for k

    std::get<0>(ijk_id) = ij_id.first;
    std::get<1>(ijk_id) = ij_id.second;
    std::get<2>(ijk_id) = cur_loc;
  }
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& POLYHEDRON
  else if (typeid(*vol_mesher) == typeid(chi_mesh::VolumeMesherExtruder))
  {
    auto extruder = (chi_mesh::VolumeMesherExtruder*)vol_mesher;

    //====================================== Create virtual cuts
    if (vol_mesher->zcuts.empty())
    {
      int num_sub_layers = extruder->vertex_layers.size()-1;

      if ((num_sub_layers%vol_mesher->options.partition_z) != 0)
      {
        chi_log.Log(LOG_ALLERROR)
          << "Number of sub-layers in extruded mesh is not divisible "
          << "by the requested number of z-partitions.";
        exit(EXIT_FAILURE);
      }

      int delta_zk = num_sub_layers/
                     vol_mesher->options.partition_z;
      for (int k=0; k<(vol_mesher->options.partition_z); k++)
      {
        int layer_index = k*delta_zk + delta_zk;
        if (layer_index > (extruder->vertex_layers.size()-1))
        {
          layer_index = (int)extruder->vertex_layers.size()-1;
          vol_mesher->zcuts.push_back(extruder->vertex_layers[layer_index]);
        }
        else
        {
          vol_mesher->zcuts.push_back(extruder->vertex_layers[layer_index]);

          if (chi_log.GetVerbosity()==LOG_0VERBOSE_2)
          {
            printf("Z-Cut %lu, %g\n",vol_mesher->zcuts.size(),
                   extruder->vertex_layers[layer_index]);
          }
        }
      }
    }


    //====================================== Scan cuts for location
    double zmin = -1.0e-16;
    double zmax =  1.0e-16;
    for (int k=0; k<(vol_mesher->zcuts.size()); k++)
    {
      zmax =  vol_mesher->zcuts[k];

      double z = cell->centroid.z;

      if (chi_log.GetVerbosity()==LOG_0VERBOSE_2)
      {
        printf("zmax = %g, zmin = %g, cell_z = %g\n",zmax,zmin,z);
      }


      if ((z > zmin) && (z < zmax))
      {
        std::get<0>(ijk_id) = ij_id.first;
        std::get<1>(ijk_id) = ij_id.second;
        std::get<2>(ijk_id) = k;

        found_partition = true;
        break;
      }
      zmin = zmax;
    }
  }//if typeid
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& POLYHEDRON
  else if (typeid(*vol_mesher) == typeid(chi_mesh::VolumeMesherPredefined3D))
  {
    if (vol_mesher->zcuts.empty())
    {
      std::get<0>(ijk_id) = ij_id.first;
      std::get<1>(ijk_id) = ij_id.second;
      std::get<2>(ijk_id) = 0;
      found_partition = true;
    }

    //====================================== Scan cuts for location
    double zmin = -1.0e16;
    double zmax =  1.0e16;
    vol_mesher->zcuts.push_back(zmax);
    for (int k=0; k<(vol_mesher->zcuts.size()); k++)
    {
      zmax =  vol_mesher->zcuts[k];

      double z = cell->centroid.z;

      if (chi_log.GetVerbosity()==LOG_0VERBOSE_2)
      {
        printf("zmax = %g, zmin = %g, cell_z = %g\n",zmax,zmin,z);
      }


      if ((z > zmin) && (z < zmax))
      {
        std::get<0>(ijk_id) = ij_id.first;
        std::get<1>(ijk_id) = ij_id.second;
        std::get<2>(ijk_id) = k;

        found_partition = true;
        break;
      }
      zmin = zmax;
    }
  }//if typeid
  else if (typeid(*vol_mesher) == typeid(chi_mesh::VolumeMesherPredefined2D))
  {
    found_partition = true;
    std::get<0>(ijk_id) = ij_id.first;
    std::get<1>(ijk_id) = ij_id.second;
    std::get<2>(ijk_id) = 0;
  }

  //================================================== Report unallocated item_id
  if (!found_partition)
  {
    chi_log.Log(LOG_ALLERROR)
    << "A cell was encountered for which "
       "no zpartition id was found";
    exit(EXIT_FAILURE);
  }

  return ijk_id;
}




//###################################################################
/**Get a list of boundary cells.*/
void chi_mesh::VolumeMesher::
 GetBoundaryCells(chi_mesh::MeshContinuum* vol_continuum)
{
  for (const auto& cell : vol_continuum->local_cells)
  {

    for (int f=0; f<cell.faces.size(); f++)
    {
      if (cell.faces[f].neighbor < 0)
      {
        vol_continuum->boundary_cell_indices.push_back(cell.global_id);
        break;
      }
    }
  }//for local item_id
  printf("Number of boundary item_id: %lu\n",vol_continuum->boundary_cell_indices.size());

  vol_continuum->ExportCellsToPython(
    "BoundaryCells.py",true,
    &vol_continuum->boundary_cell_indices);
}




//###################################################################
/**Sets material id's using a logical volume.*/
void chi_mesh::VolumeMesher::
  SetMatIDFromLogical(chi_mesh::LogicalVolume *log_vol,bool sense, int mat_id)
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Setting material id from logical volume.";
  //============================================= Get current mesh handler
  chi_mesh::MeshHandler* handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  chi_mesh::Region* cur_region = handler->region_stack.back();
  chi_mesh::MeshContinuum* vol_cont = cur_region->GetGrid();

  int num_cells_modified = 0;
  for (auto& cell : vol_cont->local_cells)
  {
    if (log_vol->Inside(cell.centroid) && sense){
      cell.material_id = mat_id;
      ++num_cells_modified;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Done setting material id from logical volume. "
    << "Number of cells modified = " << num_cells_modified << ".";
}

//###################################################################
/**Sets material id's using a logical volume.*/
void chi_mesh::VolumeMesher::
 SetBndryIDFromLogical(chi_mesh::LogicalVolume *log_vol,bool sense, int bndry_id)
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Setting boundary id from logical volume.";
  //============================================= Get current mesh handler
  chi_mesh::MeshHandler* handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  chi_mesh::Region* cur_region = handler->region_stack.back();
  chi_mesh::MeshContinuum* vol_cont = cur_region->GetGrid();

  int num_faces_modified = 0;
  for (auto& cell : vol_cont->local_cells)
  {
    for (auto& face : cell.faces)
    {
      if (face.neighbor >= 0) continue;
      if (log_vol->Inside(face.centroid) && sense){
        face.neighbor = -1*(abs(bndry_id)+1);
        ++num_faces_modified;
      }
    }
  }

  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Done setting boundary id from logical volume. "
    << "Number of faces modified = " << num_faces_modified << ".";
}