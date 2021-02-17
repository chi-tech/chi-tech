#include "diffusion_solver.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/Region/chi_region.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell.h"
#include "ChiMesh/Cell/cell_slab.h"
#include "ChiMesh/Cell/cell_polygon.h"
#include "ChiMesh/Cell/cell_polyhedron.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"
#include "ChiTimer/chi_timer.h"

#include <ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_slab.h>
#include <ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_polygon.h>
#include <ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_polyhedron.h>

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

/**Still searching for a reference for this.
 *
 * For Polygons:
 * Defined from paper  \n
 * Turcksin B, Ragusa J, "Discontinuous diffusion synthetic acceleration
 * for S_n transport on 2D arbitrary polygonal meshes", Journal of
 * Computational Physics 274, pg 356-369, 2014.\n
 * \n
 * Nv = Number of vertices. If Nv <= 4 then the perimeter parameter
 * should be replaced by edge length.*/
double chi_diffusion::Solver::HPerpendicular(chi_mesh::Cell* cell,
                                             CellPWLFEValues* fe_view,
                                             int f)
{
  double hp = 1.0;

  int Nf = cell->faces.size();
  int Nv = cell->vertex_ids.size();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
  if (cell->Type() == chi_mesh::CellType::SLAB)
  {
    auto slab_fe_view = (SlabPWLFEView*)fe_view;
    hp = slab_fe_view->h/2.0;
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  else if (cell->Type() == chi_mesh::CellType::POLYGON)
  {
    Nv = 4;
    chi_mesh::CellFace& face = cell->faces[f];

    int v0i = face.vertex_ids[0];
    int v1i = face.vertex_ids[1];

    chi_mesh::Vertex& v0 = *grid->vertices[v0i];
    chi_mesh::Vertex& v1 = *grid->vertices[v1i];

    double perimeter = (v1 - v0).Norm();

    double area  = 0.0;
    for (int i=0; i<fe_view->dofs; i++)
      area += fe_view->IntV_shapeI[i];

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
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
  else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
  {
    double volume  = 0.0;
    for (int i=0; i<fe_view->dofs; i++)
      volume += fe_view->IntV_shapeI[i];

    double area = 0.0;
    for (int fr=0; fr<Nf; fr++)
      for (int i=0; i<Nv; i++)
        area += fe_view->IntS_shapeI[i][fr];

    if (Nf == 4)                  //Tet
      hp = 3*volume/area;
    else if (Nf == 6 && Nv == 8)  //Hex
      hp = volume/area;
    else                          //Polyhedron
      hp = 6*volume/area;
  }//Polyhedron
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "Unsupported cell type in call to HPerpendicular";
    exit(EXIT_FAILURE);
  }


  return hp;
}



/**Given a global node index, returns the dof its associated on the
 * referenced cell. Polyhedron overload.*/
int chi_diffusion::Solver::MapCellDof(chi_mesh::Cell* cell, int ig)
{
  int imap = -1;
  for (int ai=0; ai < cell->vertex_ids.size(); ai++)
  {
    if (ig == cell->vertex_ids[ai])
    {
      imap = ai;
      break;
    }
  }
  return imap;
}



/**Given the face index on the current cell, finds the
 * corresponding face index on the adjacent cell.*/
int chi_diffusion::Solver::MapCellFace(chi_mesh::Cell* cur_cell,
                                       chi_mesh::Cell* adj_cell,
                                       int f)
{
  int num_face_dofs = cur_cell->faces[f].vertex_ids.size();
  int fmap = -1;
  for (int af=0; af < adj_cell->faces.size(); af++)
  {
    bool is_match = true;

    for (int fi=0; fi<num_face_dofs; fi++)
    {
      bool found = false;

      for (int afi=0; afi < adj_cell->faces[af].vertex_ids.size(); afi++)
      {
        if (cur_cell->faces[f].vertex_ids[fi] ==
            adj_cell->faces[af].vertex_ids[afi])
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
DiffusionIPCellView* chi_diffusion::Solver::GetBorderIPView(int locI,
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
CellPWLFEValues* chi_diffusion::Solver::GetBorderFEView(int locI,
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
  auto pwl_sdm = std::static_pointer_cast<SpatialDiscretization_PWL>(this->discretization);
  //============================================= Obtain the border cell info
  DiffusionIPBorderCell* cell_info =
    ip_locI_bordercell_info[locI][cell_border_index];

  //============================================= Create IP view
  DiffusionIPCellView* ip_view = new DiffusionIPCellView;
  ip_view->cell_dof_start = cell_info->cell_dof_start;
  ip_locI_borderipviews[locI][cell_border_index] = ip_view;


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
  if (cell_info->cell_type == 3)
  {
    auto cell = new chi_mesh::CellSlab;
    cell->partition_id = locI;
    cell->material_id = cell_info->cell_mat_id;

    chi_mesh::Vector3 vc;
    for (int v=0; v<cell_info->cell_dof_count; v++)
    {
      vc = vc + *grid->vertices[cell_info->v_indices[v]];
      cell->vertex_ids.push_back(cell_info->v_indices[v]);
    }
    cell->centroid = vc/cell_info->cell_dof_count;

    cell->faces.resize(cell_info->cell_face_count);
    for (int f=0; f<cell_info->cell_face_count; f++)
    {
      for (int fv=0; fv<cell_info->face_v_indices[f].size(); fv++)
        cell->faces[f].vertex_ids.push_back(cell_info->face_v_indices[f][fv]);
    }

    ip_locI_bordercells[locI][cell_border_index] = cell;

    auto fe_view = new SlabPWLFEView(cell,
                                     grid,
                                     pwl_sdm->line_quad_order_second);

    ip_locI_borderfeviews[locI][cell_border_index] = fe_view;

  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  else if (cell_info->cell_type == 4)
  {
    auto cell = new chi_mesh::CellPolygon;
    cell->partition_id = locI;
    cell->material_id = cell_info->cell_mat_id;

    chi_mesh::Vector3 vc;
    for (int v=0; v<cell_info->cell_dof_count; v++)
    {
      vc = vc + *grid->vertices[cell_info->v_indices[v]];
      cell->vertex_ids.push_back(cell_info->v_indices[v]);
    }
    cell->centroid = vc/cell_info->cell_dof_count;

    cell->faces.resize(cell_info->cell_face_count);
    for (int f=0; f<cell_info->cell_face_count; f++)
    {
      for (int fv=0; fv<cell_info->face_v_indices[f].size(); fv++)
        cell->faces[f].vertex_ids.push_back(cell_info->face_v_indices[f][fv]);
    }

    ip_locI_bordercells[locI][cell_border_index] = cell;

    auto pwl_sdm_ptr = std::static_pointer_cast<SpatialDiscretization_PWL>(discretization);


    auto fe_view =
      new PolygonPWLFEValues(cell,
                             grid,
                             pwl_sdm_ptr->tri_quad_order_second,
                             pwl_sdm_ptr->line_quad_order_second);

    fe_view->PreComputeValues();

    ip_locI_borderfeviews[locI][cell_border_index] = fe_view;

  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
  else if (cell_info->cell_type == 5)
  {
    auto cell = new chi_mesh::CellPolyhedron;
    cell->partition_id = locI;
    cell->material_id = cell_info->cell_mat_id;

    chi_mesh::Vector3 vc;
    for (int v=0; v<cell_info->cell_dof_count; v++)
    {
      vc = vc + *grid->vertices[cell_info->v_indices[v]];
      cell->vertex_ids.push_back(cell_info->v_indices[v]);
    }
    cell->centroid = vc/cell_info->cell_dof_count;

    cell->faces.resize(cell_info->cell_face_count);
    for (int f=0; f<cell_info->cell_face_count; f++)
    {
      chi_mesh::Vertex vfc;
      for (int fv=0; fv<cell_info->face_v_indices[f].size(); fv++)
      {
        cell->faces[f].vertex_ids.push_back(cell_info->face_v_indices[f][fv]);
        vfc = vfc + *grid->vertices[cell_info->face_v_indices[f][fv]];
      }
      vfc = vfc/cell_info->face_v_indices[f].size();
      cell->faces[f].centroid = vfc;

      chi_mesh::Vector3 v0fc = vfc - *grid->vertices[cell_info->face_v_indices[f][0]];
      chi_mesh::Vector3 v01 = *grid->vertices[cell_info->face_v_indices[f][1]] -
                              *grid->vertices[cell_info->face_v_indices[f][0]];

      chi_mesh::Vector3 n = v01.Cross(v0fc);
      cell->faces[f].normal = n/n.Norm();

    }//for f


    ip_locI_bordercells[locI][cell_border_index] = cell;

    auto pwl_sdm_ptr = std::static_pointer_cast<SpatialDiscretization_PWL>(discretization);

    PolyhedronPWLFEValues* fe_view =
      new PolyhedronPWLFEValues(cell,
                                grid,
                                pwl_sdm_ptr->tet_quad_order_second,
                                pwl_sdm_ptr->tri_quad_order_second);

    fe_view->PreComputeValues();

    ip_locI_borderfeviews[locI][cell_border_index] = fe_view;

  }//polyhedron
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "In function chi_diffusion::Solver::SpawnBorderCell: "
         "Unsupported cell type encountered.";
    exit(EXIT_FAILURE);
  }

}


