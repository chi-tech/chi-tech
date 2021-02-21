#include "diffusion_solver.h"


#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_slab.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_polygon.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_polyhedron.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
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
    auto v0 = *grid->vertices[cell->vertex_ids[0]];
    auto v1 = *grid->vertices[cell->vertex_ids[1]];

    hp = (v1-v0).Norm()/2.0;
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
    for (int i=0; i<fe_view->num_nodes; i++)
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
    for (int i=0; i<fe_view->num_nodes; i++)
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
double chi_diffusion::Solver::
  HPerpendicular(const chi_mesh::Cell& cell,
                 const chi_math::finite_element::UnitIntegralData& fe_intgrl_values,
                 int f)
{
  double hp = 1.0;

  int Nf = cell.faces.size();
  int Nv = cell.vertex_ids.size();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
  if (cell.Type() == chi_mesh::CellType::SLAB)
  {
    auto v0 = *grid->vertices[cell.vertex_ids[0]];
    auto v1 = *grid->vertices[cell.vertex_ids[1]];

    hp = (v1-v0).Norm()/2.0;
  }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  else if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    Nv = 4;
    const chi_mesh::CellFace& face = cell.faces[f];

    int v0i = face.vertex_ids[0];
    int v1i = face.vertex_ids[1];

    chi_mesh::Vertex& v0 = *grid->vertices[v0i];
    chi_mesh::Vertex& v1 = *grid->vertices[v1i];

    double perimeter = (v1 - v0).Norm();

    double area  = 0.0;
    for (int i=0; i<fe_intgrl_values.num_nodes; i++)
      area += fe_intgrl_values.IntV_shapeI[i];

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
  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    double volume  = 0.0;
    for (int i=0; i<fe_intgrl_values.num_nodes; i++)
      volume += fe_intgrl_values.IntV_shapeI[i];

    double area = 0.0;
    for (int fr=0; fr<Nf; fr++)
      for (int i=0; i<Nv; i++)
        area += fe_intgrl_values.IntS_shapeI[i][fr];

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



/**Given a global node index, returns the local cell-node it's associated on the
 * referenced cell. Polyhedron overload.*/
uint64_t chi_diffusion::Solver::MapCellLocalNodeIDFromGlobalID(chi_mesh::Cell* cell,
                                                               uint64_t node_global_id)
{
  size_t imap = 0;
  bool map_found = false;
  for (size_t ai=0; ai < cell->vertex_ids.size(); ai++)
  {
    if (node_global_id == cell->vertex_ids[ai])
    {
      imap = ai;
      map_found = true;
      break;
    }
  }

  if (not map_found)
    throw std::logic_error(std::string(__FUNCTION__)+": Mapping failure.");

  return imap;
}

/**Given a global node index, returns the local cell-node it's associated on the
 * referenced cell. Polyhedron overload.*/
uint64_t chi_diffusion::Solver::
  MapCellLocalNodeIDFromGlobalID(const chi_mesh::Cell& cell,
                                 uint64_t node_global_id)
{
  size_t imap = 0;
  bool map_found = false;
  for (size_t ai=0; ai < cell.vertex_ids.size(); ai++)
  {
    if (node_global_id == cell.vertex_ids[ai])
    {
      imap = ai;
      map_found = true;
      break;
    }
  }

  if (not map_found)
    throw std::logic_error(std::string(__FUNCTION__)+": Mapping failure.");

  return imap;
}



/**Given the face index on the current cell, finds the
 * corresponding face index on the adjacent cell.*/
unsigned int
chi_diffusion::Solver::MapCellFace(chi_mesh::Cell* cur_cell,
                                   chi_mesh::Cell* adj_cell,
                                   unsigned int f)
{
  const auto& ccface = cur_cell->faces[f]; //current cell face
  std::set<uint64_t> ccface_vids;
  for (auto vid : ccface.vertex_ids) ccface_vids.insert(vid);

  size_t fmap;
  bool map_found = false;
  for (size_t af=0; af < adj_cell->faces.size(); af++)
  {
    const auto& acface = adj_cell->faces[af]; //adjacent cell face

    std::set<uint64_t> acface_vids;
    for (auto vid : acface.vertex_ids) acface_vids.insert(vid);

    if (acface_vids == ccface_vids)
    {
      fmap = af;
      map_found = true;
      break;
    }
  }//for adj faces

  if (not map_found)
    throw std::logic_error(std::string(__FUNCTION__)+": Mapping failure.");

  return (unsigned int)fmap;
}

/**Given the face index on the current cell, finds the
 * corresponding face index on the adjacent cell.*/
unsigned int
chi_diffusion::Solver::MapCellFace(const chi_mesh::Cell& cur_cell,
                                   const chi_mesh::Cell& adj_cell,
                                   unsigned int f)
{
  const auto& ccface = cur_cell.faces[f]; //current cell face
  std::set<uint64_t> ccface_vids;
  for (auto vid : ccface.vertex_ids) ccface_vids.insert(vid);

  size_t fmap;
  bool map_found = false;
  for (size_t af=0; af < adj_cell.faces.size(); af++)
  {
    const auto& acface = adj_cell.faces[af]; //adjacent cell face

    std::set<uint64_t> acface_vids;
    for (auto vid : acface.vertex_ids) acface_vids.insert(vid);

    if (acface_vids == ccface_vids)
    {
      fmap = af;
      map_found = true;
      break;
    }
  }//for adj faces

  if (not map_found)
    throw std::logic_error(std::string(__FUNCTION__)+": Mapping failure.");

  return (unsigned int)fmap;
}


