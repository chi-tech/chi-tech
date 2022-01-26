#include "chi_volumemesher.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Region/chi_region.h"
#include "ChiMesh/SurfaceMesh/chi_surfacemesh.h"
#include "ChiMesh/SurfaceMesher/Predefined/surfmesher_predefined.h"
#include "ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h"
#include "ChiMesh/VolumeMesher/PredefinedUnpartitioned/volmesher_predefunpart.h"
#include "ChiMesh/Cell/cell.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiMPI/chi_mpi.h"
extern ChiMPI& chi_mpi;

#include <ChiTimer/chi_timer.h>
extern ChiTimer chi_program_timer;

//###################################################################
/**Creates 2D polygon cells for each face of a surface mesh.*/
void chi_mesh::VolumeMesher::
  AddContinuumToRegion(chi_mesh::MeshContinuumPtr& grid,
                       chi_mesh::Region& region)
{
  region.volume_mesh_continua.push_back(grid);
}

//###################################################################
/**Creates 2D polygon cells for each face of a surface mesh.*/
void chi_mesh::VolumeMesher::
  CreatePolygonCells(chi_mesh::SurfaceMesh *surface_mesh,
                     chi_mesh::MeshContinuumPtr& vol_continuum,
                     bool delete_surface_mesh_elements,
                     bool force_local)
{
  //============================================= Get current mesh handler
  chi_mesh::MeshHandler* handler = chi_mesh::GetCurrentHandler();

  //============================================= Copy nodes
  {
    uint64_t id = 0;
    for (auto& vertex : surface_mesh->vertices)
      vol_continuum->vertices.Insert(id++, vertex);
  }

  //============================================= Delete nodes
  if (delete_surface_mesh_elements)
    surface_mesh->vertices = std::vector<chi_mesh::Vertex>(0);

  //================================== Check if already have material ids
  bool have_mat_ids = false;
  if (surface_mesh->physical_region_map.size() == surface_mesh->poly_faces.size())
    have_mat_ids = true;
  else if (surface_mesh->physical_region_map.size())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Material id's specified, but are inconsistent with mesh.";
    exit(EXIT_FAILURE);
  }


  //============================================= Process faces
  unsigned int num_cells = 0;
  for (auto& face : surface_mesh->faces)
  {
    auto cell = new chi_mesh::Cell(CellType::POLYGON,CellType::TRIANGLE);

    for (int k=0;k<3;k++)
    {
      cell->vertex_ids.push_back(face.v_index[k]);

      chi_mesh::CellFace new_face;

      new_face.vertex_ids.push_back(face.e_index[k][0]);
      new_face.vertex_ids.push_back(face.e_index[k][1]);


      const auto& v0 = vol_continuum->vertices[face.e_index[k][0]];
      const auto& v1 = vol_continuum->vertices[face.e_index[k][1]];
      new_face.centroid = v0*0.5 + v1*0.5;

      chi_mesh::Vector3 vk = chi_mesh::Vector3(0.0, 0.0, 1.0);

      chi_mesh::Vector3 va = v1 - v0;
      chi_mesh::Vector3 vn = va.Cross(vk);
      vn = vn/vn.Norm();
      new_face.normal = vn;

      if (face.e_index[k][2]>=0)
      {
        new_face.neighbor_id = face.e_index[k][2];
        new_face.has_neighbor = true;
      }

      cell->faces.push_back(new_face);

      cell->centroid = cell->centroid + surface_mesh->vertices[face.v_index[k]];
    }
    cell->centroid = cell->centroid/3;

    //====================================== Compute xy partition id
    auto xy_partition_indices = GetCellXYPartitionID(cell);
    cell->partition_id = xy_partition_indices.second*
                         handler->volume_mesher->options.partition_x +
                         xy_partition_indices.first;

    if (force_local)
      cell->partition_id = chi_mpi.location_id;

    cell->global_id = num_cells;

    vol_continuum->cells.push_back(cell); ++num_cells;
  }

  for (auto face : surface_mesh->poly_faces)
  {
    CellType sub_type = CellType::POLYGON;

    const size_t num_verts = face->v_indices.size();
    if      (num_verts == 3) sub_type = CellType::TRIANGLE;
    else if (num_verts == 4) sub_type = CellType::QUADRILATERAL;

    auto cell = new chi_mesh::Cell(CellType::POLYGON, sub_type);

    //====================================== Copy vertices
    for (auto vid : face->v_indices)
    {
      cell->vertex_ids.push_back(vid);
      cell->centroid = cell->centroid + vol_continuum->vertices[vid];
    }
    cell->centroid = cell->centroid/cell->vertex_ids.size();

    //====================================== Copy edges
    for (auto src_side : face->edges)
    {
      chi_mesh::CellFace new_face;

      new_face.vertex_ids.push_back(src_side[0]);
      new_face.vertex_ids.push_back(src_side[1]);

      const auto& v0 = vol_continuum->vertices[src_side[0]];
      const auto& v1 = vol_continuum->vertices[src_side[1]];
      new_face.centroid = v0*0.5 + v1*0.5;
      chi_mesh::Vector3 vk = chi_mesh::Vector3(0.0, 0.0, 1.0);

      chi_mesh::Vector3 va = v1 - v0;
      chi_mesh::Vector3 vn = va.Cross(vk);
      vn = vn/vn.Norm();
      new_face.normal = vn;

      if (src_side[2] >= 0)
      {
        new_face.neighbor_id = src_side[2];
        new_face.has_neighbor = true;
      }
      else
        new_face.neighbor_id = 0;

      cell->faces.push_back(new_face);
    }

    //====================================== Compute partition id
    auto xy_partition_indices = GetCellXYPartitionID(cell);
    cell->partition_id = xy_partition_indices.second*
                         handler->volume_mesher->options.partition_x +
                         xy_partition_indices.first;

    if (force_local)
      cell->partition_id = chi_mpi.location_id;

    cell->global_id = num_cells;

    if (have_mat_ids)
      cell->material_id = surface_mesh->physical_region_map[num_cells];

    vol_continuum->cells.push_back(cell);
    ++num_cells;

    if (delete_surface_mesh_elements)
      delete face;
  }

  if (delete_surface_mesh_elements)
    surface_mesh->poly_faces.clear();

}

//###################################################################
/**Creates 2D polygon cells for each face of a surface mesh.*/
void chi_mesh::VolumeMesher::
  CreatePolygonCells(const chi_mesh::UnpartitionedMesh& umesh,
                     chi_mesh::MeshContinuumPtr& grid)
{
  //=================================== Copy nodes
  {
    uint64_t id = 0;
    for (const auto& vertex : umesh.vertices)
      grid->vertices.Insert(id++, vertex);
  }

  size_t num_cells=0;
  for (auto& raw_cell : umesh.raw_cells)
  {
    // Check valid template cell
    if (raw_cell->type != chi_mesh::CellType::POLYGON)
    {
      chi_log.Log(LOG_ALLERROR)
        << "chi_mesh::VolumeMesher::CreatePolygonCells "
           "called with a cell not being of primary type"
           " chi_mesh::CellType::POLYGON.";
      exit(EXIT_FAILURE);
    }

    //====================================== Make cell
    auto cell = new chi_mesh::Cell(CellType::POLYGON, raw_cell->sub_type);

    cell->global_id = num_cells;
    cell->local_id  = num_cells;
    cell->partition_id = chi_mpi.location_id;

    cell->centroid    = raw_cell->centroid;
    cell->material_id = raw_cell->material_id;
    cell->vertex_ids  = raw_cell->vertex_ids;

    // Copy faces + compute face centroid and normal
    const chi_mesh::Vector3 khat(0.0, 0.0, 1.0);
    for (auto& raw_face : raw_cell->faces)
    {
      chi_mesh::CellFace new_face;

      new_face.vertex_ids  = raw_face.vertex_ids;

      const auto& v0 = grid->vertices[new_face.vertex_ids[0]];
      const auto& v1 = grid->vertices[new_face.vertex_ids[1]];
      new_face.centroid = v0*0.5 + v1*0.5;

      chi_mesh::Vector3 va = v1 - v0;
      chi_mesh::Vector3 vn = va.Cross(khat);
      vn = vn/vn.Norm();
      new_face.normal = vn;

      new_face.has_neighbor = raw_face.has_neighbor;
      new_face.neighbor_id = raw_face.neighbor;

      cell->faces.push_back(new_face);
    }

    //====================================== Push to grid
    grid->cells.push_back(cell);
    ++num_cells;
  }//for raw_cell
}

//###################################################################
/**Obtains the xy partition IDs of a cell.
 * Cell xy_partition ids are obtained from
 * the surface mesher.*/
std::pair<int,int> chi_mesh::VolumeMesher::
 GetCellXYPartitionID(chi_mesh::Cell *cell)
{
  std::pair<int,int> ij_id(0,0);

  if (chi_mpi.process_count == 1){return ij_id;}

  //================================================== Get the current handler
  auto mesh_handler = chi_mesh::GetCurrentHandler();
  auto vol_mesher = mesh_handler->volume_mesher;

//====================================== Sanity check on partitioning
  size_t num_x_subsets = vol_mesher->options.xcuts.size()+1;
  size_t num_y_subsets = vol_mesher->options.ycuts.size()+1;

  size_t x_remainder = num_x_subsets%vol_mesher->options.partition_x;
  size_t y_remainder = num_y_subsets%vol_mesher->options.partition_y;

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

  size_t subsets_per_partitionx = num_x_subsets/vol_mesher->options.partition_x;
  size_t subsets_per_partitiony = num_y_subsets/vol_mesher->options.partition_y;

  //====================================== Determine x-partition
  int x=-1;
  int xcount=-1;
  for (size_t i =  subsets_per_partitionx-1;
       i <  vol_mesher->options.xcuts.size();
       i += subsets_per_partitionx)
  {
    xcount++;
    if (cell->centroid.x <= vol_mesher->options.xcuts[i])
    {
      x = xcount;
      break;
    }
  }
  if (x<0)
  {
    x = vol_mesher->options.partition_x-1;
  }

  //====================================== Determine y-partition
  int y=-1;
  int ycount=-1;
  for (size_t i =  subsets_per_partitiony-1;
       i <  vol_mesher->options.ycuts.size();
       i += subsets_per_partitiony)
  {
    ycount++;
    if (cell->centroid.y <= vol_mesher->options.ycuts[i])
    {
      y = ycount;
      break;
    }
  }
  if (y<0)
  {
    y = vol_mesher->options.partition_y - 1;
  }

  //====================================== Set partitioning
  ij_id.first = x;
  ij_id.second= y;

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

  if (vol_mesher->options.partition_z == 1)
  {
    found_partition = true;
    std::get<0>(ijk_id) = ij_id.first;
    std::get<1>(ijk_id) = ij_id.second;
    std::get<2>(ijk_id) = 0;
  }
  else if (typeid(*vol_mesher) == typeid(chi_mesh::VolumeMesherExtruder))
  {
    auto extruder = (chi_mesh::VolumeMesherExtruder*)vol_mesher;

    //====================================== Create virtual cuts
    if (vol_mesher->options.zcuts.empty())
    {
      size_t num_sub_layers = extruder->vertex_layers.size()-1;

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
          vol_mesher->options.zcuts.push_back(extruder->vertex_layers[layer_index]);
        }
        else
        {
          vol_mesher->options.zcuts.push_back(extruder->vertex_layers[layer_index]);

          if (chi_log.GetVerbosity()==LOG_0VERBOSE_2)
          {
            printf("Z-Cut %lu, %g\n",vol_mesher->options.zcuts.size(),
                   extruder->vertex_layers[layer_index]);
          }
        }
      }
    }


    //====================================== Scan cuts for location
    double zmin = -1.0e-16;
    for (int k=0; k<(vol_mesher->options.zcuts.size()); k++)
    {
      double zmax =  vol_mesher->options.zcuts[k];

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
  else if (typeid(*vol_mesher) == typeid(chi_mesh::VolumeMesherPredefinedUnpartitioned))
  {
    if (vol_mesher->options.zcuts.empty())
    {
      throw std::invalid_argument("Cell z-partitioning cannot be determined "
                                  "because no z-cuts are supplied to volume "
                                  "mesher.");
    }

    //====================================== Scan cuts for location
    std::vector<double> temp_zcuts = vol_mesher->options.zcuts;
    double zmin = -1.0e16;
    double zmax =  1.0e16;
    temp_zcuts.push_back(zmax);
    for (int k=0; k<(temp_zcuts.size()); k++)
    {
      zmax =  temp_zcuts[k];

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
    }//for k
  }//if typeid

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
  chi_mesh::MeshContinuumPtr vol_cont = cur_region->GetGrid();

  int num_cells_modified = 0;
  for (auto& cell : vol_cont->local_cells)
  {
    if (log_vol->Inside(cell.centroid) && sense){
      cell.material_id = mat_id;
      ++num_cells_modified;
    }
  }

  const auto& ghost_ids = vol_cont->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = vol_cont->cells[ghost_id];
    if (log_vol->Inside(cell.centroid) && sense)
      cell.material_id = mat_id;
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
  chi_mesh::MeshContinuumPtr vol_cont = cur_region->GetGrid();

  int num_faces_modified = 0;
  for (auto& cell : vol_cont->local_cells)
  {
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor) continue;
      if (log_vol->Inside(face.centroid) && sense){
        face.neighbor_id = abs(bndry_id);
        ++num_faces_modified;
      }
    }
  }

  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Done setting boundary id from logical volume. "
    << "Number of faces modified = " << num_faces_modified << ".";
}

//###################################################################
/**Sets material id's for all cells to the specified material id.*/
void chi_mesh::VolumeMesher::
  SetMatIDToAll(int mat_id)
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Setting material id " << mat_id << "to all cells.";

  //============================================= Get current mesh handler
  auto handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  auto cur_region = handler->region_stack.back();
  auto vol_cont = cur_region->GetGrid();

  for (auto& cell : vol_cont->local_cells)
    cell.material_id = mat_id;

  const auto& ghost_ids = vol_cont->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
    vol_cont->cells[ghost_id].material_id = mat_id;

  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Done setting material id " << mat_id << " to all cells";
}

//###################################################################
/**Sets boundary numbers on boundaries orthogonal to the cardinal directions
 * as xmax=0, xmin=1, ymax=2, ymin=3, zmax=4, zmin=5.*/
void chi_mesh::VolumeMesher::
  SetupOrthogonalBoundaries()
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Setting orthogonal boundaries.";

  //============================================= Get current mesh handler
  auto handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  auto cur_region = handler->region_stack.back();
  auto vol_cont = cur_region->GetGrid();

  const chi_mesh::Vector3 ihat(1.0, 0.0, 0.0);
  const chi_mesh::Vector3 jhat(0.0, 1.0, 0.0);
  const chi_mesh::Vector3 khat(0.0, 0.0, 1.0);

  for (auto& cell : vol_cont->local_cells)
    for (auto& face : cell.faces)
      if (not face.has_neighbor)
      {
        chi_mesh::Vector3& n = face.normal;

        int boundary_id = -1;
        if      (n.Dot(ihat)>0.999)  boundary_id = 0;
        else if (n.Dot(ihat)<-0.999) boundary_id = 1;
        else if (n.Dot(jhat)> 0.999) boundary_id = 2;
        else if (n.Dot(jhat)<-0.999) boundary_id = 3;
        else if (n.Dot(khat)> 0.999) boundary_id = 4;
        else if (n.Dot(khat)<-0.999) boundary_id = 5;

        if (boundary_id >= 0) face.neighbor_id = boundary_id;
      }//if bndry

  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Done setting orthogonal boundaries.";
}
