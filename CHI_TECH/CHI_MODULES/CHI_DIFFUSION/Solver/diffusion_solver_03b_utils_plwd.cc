#include "diffusion_solver.h"
#include "../../../CHI_MESH/CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../../CHI_MESH/CHI_REGION/chi_region.h"
#include "../../../CHI_MESH/CHI_MESHCONTINUUM/chi_meshcontinuum.h"
#include "../../../CHI_MESH/CHI_CELL/cell.h"
#include "../../../CHI_MESH/CHI_CELL/cell_slab.h"
#include "../../../CHI_MESH/CHI_CELL/cell_polygon.h"
#include "../../../CHI_MESH/CHI_CELL/cell_polyhedron.h"
#include "../../../CHI_MESH/CHI_VOLUMEMESHER/chi_volumemesher.h"
#include "../../../CHI_TIMER/chi_timer.h"

#include <ChiMath/SpatialDiscretization/PiecewiseLinear/CellViews/pwl_slab.h>
#include <ChiMath/SpatialDiscretization/PiecewiseLinear/CellViews/pwl_polygon.h>
#include <ChiMath/SpatialDiscretization/PiecewiseLinear/CellViews/pwl_polyhedron.h>

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <chi_log.h>
#include <chi_mpi.h>

extern CHI_LOG chi_log;
extern CHI_MPI chi_mpi;



//###################################################################
/**Reorders nodes for better parrallel communication during matrix
 * assembly specific to PWLD methods.*/
void chi_diffusion::Solver::ReorderNodesPWLD()
{
  CHI_TIMER t_stage[6];

  t_stage[0].Reset();
  //================================================== Get reference to continuum
  auto handler = chi_mesh::GetCurrentHandler();
  auto region  = handler->region_stack.back();
  auto vol_continuum = region->volume_mesh_continua.back();

  auto mesher = handler->volume_mesher;

  //================================================== Get local DOF count
  pwld_local_dof_count=0;
  size_t num_loc_cells = vol_continuum->local_cell_glob_indices.size();
  for (int lc=0; lc<num_loc_cells; lc++)
  {
    int cell_glob_index = vol_continuum->local_cell_glob_indices[lc];
    auto cell = vol_continuum->cells[cell_glob_index];

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      pwld_local_dof_count += 2;
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      auto poly_cell = (chi_mesh::CellPolygon*)cell;
      pwld_local_dof_count += poly_cell->v_indices.size();
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
      pwld_local_dof_count += polyh_cell->v_indices.size();
    }
  }

  //================================================== Get global DOF count
  pwld_global_dof_count=0;
  MPI_Allreduce(&pwld_local_dof_count,    //Send buffer
                &pwld_global_dof_count,   //Recv buffer
                1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  //================================================== Ring communicate DOF start
  pwld_local_dof_start = 0;
  if (chi_mpi.location_id != 0)
  {
    MPI_Recv(&pwld_local_dof_start,
             1,MPI_INT,              //Count and type
             chi_mpi.location_id-1,  //Source
             111,                    //Tag
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }

  if (chi_mpi.location_id != (chi_mpi.process_count-1))
  {
    int next_loc_start = pwld_local_dof_start+pwld_local_dof_count;
    MPI_Send(&next_loc_start,
             1,MPI_INT,
             chi_mpi.location_id+1,
             111,
             MPI_COMM_WORLD);
  }

  chi_log.Log(LOG_ALLVERBOSE_2)
   << "Local dof count, start, total "
   << pwld_local_dof_count << " "
   << pwld_local_dof_start << " "
   << pwld_global_dof_count;

}




//###################################################################
/**Maps a border cell. Returns the matrix row location for the given
 * vertex index.*/
int chi_diffusion::Solver::MapBorderCell(int locI, int neighbor, int vglob_i)
{
  for (int c=0; c<ip_locI_bordercell_info[locI].size(); c++)
  {
    DIFFUSION_IP_BORDERCELL* border_cell = ip_locI_bordercell_info[locI][c];

    if (border_cell->cell_glob_index == neighbor)
    {
      int ref_dof=-1;
      for (int v=0; v<border_cell->v_indices.size(); v++)
      {
        if (border_cell->v_indices[v] == vglob_i)
        {
          return border_cell->cell_dof_start + v;
        }
      }

      chi_log.Log(LOG_ALLERROR)
        << "Failed call to chi_diffusion::Solver::MapBorderCell. "
        << "Failed to find associated vertex.";
      exit(EXIT_FAILURE);
    }
  }

  chi_log.Log(LOG_ALLERROR)
    << "Failed call to chi_diffusion::Solver::MapBorderCell. "
    << "Failed to find cell associated with neighbor.";
  exit(EXIT_FAILURE);

  return -1;
}


/**Defined from paper  \n
 * Turcksin B, Ragusa J, "Discontinuous diffusion synthetic acceleration
 * for S_n transport on 2D arbitrary polygonal meshes", Journal of
 * Computational Physics 274, pg 356-369, 2014.\n
 * \n
 * Nv = Number of vertices. If Nv <= 4 then the perimeter parameter
 * should be replaced by edge length.*/
double chi_diffusion::Solver::HPerpendicularPoly(int Nv, double area, double perimeter)
{
  double hp = 0.0;

  if (Nv == 3)
    hp = 2*area/perimeter;
  else if (Nv == 4)
    hp = area/perimeter;
  else //Nv > 4
  {
    if (Nv%2 == 0)
      hp = 4*area/perimeter;
    else
    {
      hp = 2*area/perimeter;
      hp += sqrt(2*area / Nv*sin(2*M_PI/Nv));
    }
  }

  return hp;
}

/**Still searching for a reference for this.*/
double chi_diffusion::Solver::HPerpendicularPolyH(int Nf, int Nv, double volume, double area)
{
  double hp = 0.0;

  if (Nf == 4)                  //Tet
    hp = 3*volume/area;
  else if (Nf == 6 && Nv == 8)  //Hex
    hp = volume/area;
  else                          //Polyhedron
    hp = 6*volume/area;

  return hp;
}

/**Given a global node index, returns the dof its associated on the
 * referenced cell. Slab overload.*/
int chi_diffusion::Solver::MapCellDof(chi_mesh::CellSlab* slab_cell, int ig)
{
  int imap = -1;
  for (int ai=0; ai<2; ai++)
  {
    if (ig == slab_cell->v_indices[ai])
    {
      imap = ai;
      break;
    }
  }
  return imap;
}

/**Given a global node index, returns the dof its associated on the
 * referenced cell. Polygon overload.*/
int chi_diffusion::Solver::MapCellDof(chi_mesh::CellPolygon* poly_cell, int ig)
{
  int imap = -1;
  for (int ai=0; ai<poly_cell->v_indices.size(); ai++)
  {
    if (ig == poly_cell->v_indices[ai])
    {
      imap = ai;
      break;
    }
  }
  return imap;
}

/**Given a global node index, returns the dof its associated on the
 * referenced cell. Polyhedron overload.*/
int chi_diffusion::Solver::MapCellDof(chi_mesh::CellPolyhedron* polyh_cell, int ig)
{
  int imap = -1;
  for (int ai=0; ai<polyh_cell->v_indices.size(); ai++)
  {
    if (ig == polyh_cell->v_indices[ai])
    {
      imap = ai;
      break;
    }
  }
  return imap;
}

/**Given the face index on the current cell, finds the
 * corresponding face index on the adjacent cell.*/
int chi_diffusion::Solver::MapCellFace(chi_mesh::CellPolyhedron* polyh_cell,
                                       chi_mesh::CellPolyhedron* adjph_cell,
                                       int f)
{
  int num_face_dofs = polyh_cell->faces[f]->v_indices.size();
  int fmap = -1;
  for (int af=0; af<adjph_cell->faces.size(); af++)
  {
    bool is_match = true;

    for (int fi=0; fi<num_face_dofs; fi++)
    {
      bool found = false;

      for (int afi=0; afi<adjph_cell->faces[af]->v_indices.size(); afi++)
      {
        if (polyh_cell->faces[f]->v_indices[fi] ==
          adjph_cell->faces[af]->v_indices[afi])
        {
          found = true;
          break;
        }
      }//for adj face verts

      if (!found)
      {
        is_match = false;
        break;
      }
    }//for cur face verts

    if (is_match)
      fmap = af;
  }//for adj faces

  if (fmap<0)
  {
    chi_log.Log(LOG_ALL)
      << "Error fmap";
    exit(EXIT_FAILURE);
  }

  return fmap;
}




/**Obtains a reference to a Interior Penalty view of a cell.*/
DIFFUSION_IP_VIEW* chi_diffusion::Solver::GetBorderIPView(int locI,
                                                          int cell_glob_index)
{
  int cell_border_index=-1;
  for (int c=0; c<ip_locI_bordercell_info[locI].size(); c++)
  {
    if (ip_locI_bordercell_info[locI][c]->cell_glob_index == cell_glob_index)
      cell_border_index = c;
  }

  if (cell_border_index<0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "In function chi_diffusion::Solver::GetBorderIPView: "
         "Cell border index mapping failed.";
    exit(EXIT_FAILURE);
  }

  if (ip_locI_borderipviews[locI][cell_border_index] == nullptr)
  {
    SpawnBorderCell(locI,cell_border_index);
    return ip_locI_borderipviews[locI][cell_border_index];
  }
  else
    return ip_locI_borderipviews[locI][cell_border_index];
}

/**Obtains a reference to a Finite Element view of a cell.*/
CellFEView* chi_diffusion::Solver::GetBorderFEView(int locI,
                                                   int cell_glob_index)
{
  int cell_border_index=-1;
  for (int c=0; c<ip_locI_bordercell_info[locI].size(); c++)
  {
    if (ip_locI_bordercell_info[locI][c]->cell_glob_index == cell_glob_index)
      cell_border_index = c;
  }

  if (cell_border_index<0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "In function chi_diffusion::Solver::GetBorderFEView: "
         "Cell border index mapping failed.";
    exit(EXIT_FAILURE);
  }

  if (ip_locI_borderipviews[locI][cell_border_index] == nullptr)
  {
    SpawnBorderCell(locI,cell_border_index);
    return ip_locI_borderfeviews[locI][cell_border_index];
  }
  else
    return ip_locI_borderfeviews[locI][cell_border_index];
}

/**Obtains a reference to a cell.*/
chi_mesh::Cell* chi_diffusion::Solver::GetBorderCell(int locI,
                                                     int cell_glob_index)
{
  int cell_border_index=-1;
  for (int c=0; c<ip_locI_bordercell_info[locI].size(); c++)
  {
    if (ip_locI_bordercell_info[locI][c]->cell_glob_index == cell_glob_index)
      cell_border_index = c;
  }

  if (cell_border_index<0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "In function chi_diffusion::Solver::GetBorderCell: "
         "Cell border index mapping failed.";
    exit(EXIT_FAILURE);
  }

  if (ip_locI_borderipviews[locI][cell_border_index] == nullptr)
  {
    SpawnBorderCell(locI,cell_border_index);
    return ip_locI_bordercells[locI][cell_border_index];
  }
  else
    return ip_locI_bordercells[locI][cell_border_index];
}






/**This function can be called from a number of places and is used to
 * generate the necessary information for the successful implementation
 * of MIP.*/
void chi_diffusion::Solver::SpawnBorderCell(int locI, int cell_border_index)
{
  //============================================= Obtain the border cell info
  DIFFUSION_IP_BORDERCELL* cell_info =
    ip_locI_bordercell_info[locI][cell_border_index];

  //============================================= Create IP view
  DIFFUSION_IP_VIEW* ip_view = new DIFFUSION_IP_VIEW;
  ip_view->cell_dof_start = cell_info->cell_dof_start;
  ip_locI_borderipviews[locI][cell_border_index] = ip_view;

  //============================================= Create cell
  if (cell_info->cell_type == 0)
  {
    chi_mesh::CellSlab* cell = new chi_mesh::CellSlab;
    cell->partition_id = locI;
    cell->material_id = cell_info->cell_mat_id;

    chi_mesh::Vector vc;
    for (int v=0; v<cell_info->cell_dof_count; v++)
    {
      vc = vc + *grid->nodes[cell_info->v_indices[v]];
      cell->v_indices[v] = cell_info->v_indices[v];
    }
    cell->centroid = vc/cell_info->cell_dof_count;

    ip_locI_bordercells[locI][cell_border_index] = cell;

    SlabFEView* fe_view = new SlabFEView(cell, grid);

    ip_locI_borderfeviews[locI][cell_border_index] = fe_view;
  }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  else if (cell_info->cell_type == 1)
  {
    chi_mesh::CellPolygon* cell = new chi_mesh::CellPolygon;
    cell->partition_id = locI;
    cell->material_id = cell_info->cell_mat_id;

    chi_mesh::Vector vc;
    for (int v=0; v<cell_info->cell_dof_count; v++)
    {
      vc = vc + *grid->nodes[cell_info->v_indices[v]];
      cell->v_indices.push_back(cell_info->v_indices[v]);
    }
    cell->centroid = vc/cell_info->cell_dof_count;

    cell->edges.resize(cell_info->cell_face_count);
    for (int f=0; f<cell_info->cell_face_count; f++)
    {
      cell->edges[f] = new int[2];
      for (int fv=0; fv<cell_info->face_v_indices[f].size(); fv++)
        cell->edges[f][fv] = cell_info->face_v_indices[f][fv];
    }

    ip_locI_bordercells[locI][cell_border_index] = cell;

    PolygonFEView* fe_view =
      new PolygonFEView(cell, grid, (SpatialDiscretization_PWL*)discretization);

    fe_view->PreCompute();

    ip_locI_borderfeviews[locI][cell_border_index] = fe_view;

  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
  else if (cell_info->cell_type == 2)
  {
    chi_mesh::CellPolyhedron* cell = new chi_mesh::CellPolyhedron;
    cell->partition_id = locI;
    cell->material_id = cell_info->cell_mat_id;

    chi_mesh::Vector vc;
    for (int v=0; v<cell_info->cell_dof_count; v++)
    {
      vc = vc + *grid->nodes[cell_info->v_indices[v]];
      cell->v_indices.push_back(cell_info->v_indices[v]);
    }
    cell->centroid = vc/cell_info->cell_dof_count;

    cell->faces.resize(cell_info->cell_face_count);
    for (int f=0; f<cell_info->cell_face_count; f++)
    {
      cell->faces[f] = new chi_mesh::PolyFace;
      chi_mesh::Vertex vfc;
      for (int fv=0; fv<cell_info->face_v_indices[f].size(); fv++)
      {
        cell->faces[f]->v_indices.push_back(cell_info->face_v_indices[f][fv]);
        vfc = vfc + *grid->nodes[cell_info->face_v_indices[f][fv]];
      }
      vfc = vfc/cell_info->face_v_indices[f].size();
      cell->faces[f]->face_centroid = vfc;

      chi_mesh::Vector v0fc = vfc - *grid->nodes[cell_info->face_v_indices[f][0]];
      chi_mesh::Vector v01 = *grid->nodes[cell_info->face_v_indices[f][1]] -
                             *grid->nodes[cell_info->face_v_indices[f][0]];

      chi_mesh::Vector n = v01.Cross(v0fc);
      cell->faces[f]->geometric_normal = n/n.Norm();

      //Populate edges
      int edge_count = cell->faces[f]->v_indices.size();
      cell->faces[f]->edges.resize(edge_count);
      for (int e=0; e<edge_count; e++)
      {
        cell->faces[f]->edges[e] = new int[4];
        if (e != (edge_count-1))
        {
          cell->faces[f]->edges[e][0] = cell->faces[f]->v_indices[e];
          cell->faces[f]->edges[e][1] = cell->faces[f]->v_indices[e+1];
        } else
        {
          cell->faces[f]->edges[e][0] = cell->faces[f]->v_indices[e];
          cell->faces[f]->edges[e][1] = cell->faces[f]->v_indices[0];
        }
      }

    }//for f


    ip_locI_bordercells[locI][cell_border_index] = cell;

    PolyhedronFEView* fe_view =
      new PolyhedronFEView(cell,grid,(SpatialDiscretization_PWL*)discretization);

    fe_view->PreCompute();

    ip_locI_borderfeviews[locI][cell_border_index] = fe_view;

  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "In function chi_diffusion::Solver::SpawnBorderCell: "
         "Unsupported cell type encountered.";
    exit(EXIT_FAILURE);
  }

}


