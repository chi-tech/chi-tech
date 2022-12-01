#include "diffusion_mip.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"

#include "ChiMesh/chi_mesh.h"

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
double lbs::acceleration::DiffusionMIPSolver::
  HPerpendicular(const chi_mesh::Cell& cell,
                 unsigned int f)
{
  const auto& cell_mapping = m_sdm.GetCellMapping(cell);
  double hp;

  size_t Nf = cell.faces.size();
  size_t Nv = cell.vertex_ids.size();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
  if (cell.Type() == chi_mesh::CellType::SLAB)
    hp = cell_mapping.CellVolume()/2.0;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  else if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    const chi_mesh::CellFace& face = cell.faces[f];

    uint64_t v0i = face.vertex_ids[0];
    uint64_t v1i = face.vertex_ids[1];

    const auto& v0 = m_grid.vertices[v0i];
    const auto& v1 = m_grid.vertices[v1i];

    double perimeter = (v1 - v0).Norm();

    double area  = cell_mapping.CellVolume();

    hp = area/perimeter;

//    if (Nv == 3)
//      hp = 2*area/perimeter;
//    else if (Nv == 4)
//      hp = area/perimeter;
//    else //Nv > 4
//    {
//      if (Nv%2 == 0)
//        hp = 4*area/perimeter;
//      else
//      {
//        hp = 2*area/perimeter;
//        hp += sqrt(2*area / Nv*sin(2*M_PI/Nv));
//      }
//    }
  }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    double volume  = cell_mapping.CellVolume();

    double area = 0.0;
    for (size_t fr=0; fr<Nf; ++fr)
      area += cell_mapping.FaceArea(fr);

    if (Nf == 4)                  //Tet
      hp = 3*volume/area;
    else if (Nf == 6 && Nv == 8)  //Hex
      hp = volume/area;
    else                          //Polyhedron
      hp = 6*volume/area;
  }//Polyhedron
  else
    throw std::logic_error(
      "lbs::acceleration::DiffusionMIPSolver::HPerpendicular: "
      "Unsupported cell type in call to HPerpendicular");

  return hp;
}

int lbs::acceleration::DiffusionMIPSolver::
  MapFaceNodeDisc(const chi_mesh::Cell& cur_cell,
                  const chi_mesh::Cell& adj_cell,
                  const size_t f,
                  const size_t fi,
                  const double epsilon/*=1.0e-12*/)
{
  const auto& cur_cell_mapping = m_sdm.GetCellMapping(cur_cell);
  const auto& adj_cell_mapping = m_sdm.GetCellMapping(adj_cell);

  const auto cur_cell_node_locs = cur_cell_mapping.GetNodeLocations();
  const auto adj_cell_node_locs = adj_cell_mapping.GetNodeLocations();

  const int i = cur_cell_mapping.MapFaceNode(f, fi);
  const auto& node_i_loc = cur_cell_node_locs[i];

  const size_t af = chi_mesh::MeshContinuum::MapCellFace(cur_cell,adj_cell,f);
  const size_t adj_face_num_nodes = adj_cell_mapping.NumFaceNodes(af);

  for (size_t fj=0; fj<adj_face_num_nodes; ++fj)
  {
    const int j = adj_cell_mapping.MapFaceNode(af,fj);
    if ((node_i_loc - adj_cell_node_locs[j]).NormSquare() < epsilon)
      return j;
  }

  throw std::logic_error(
    "lbs::acceleration::DiffusionMIPSolver::MapFaceNodeDisc: Mapping failure.");
}

int lbs::acceleration::DiffusionMIPSolver::
  MapFaceNodeDisc(const chi_mesh::Cell& cur_cell,
                  const chi_mesh::Cell& adj_cell,
                  const std::vector<chi_mesh::Vector3>& cc_node_locs,
                  const std::vector<chi_mesh::Vector3>& ac_node_locs,
                  size_t ccf, size_t acf,
                  size_t ccfi,
                  double epsilon/*=1.0e-12*/)
{
  const auto& cur_cell_mapping = m_sdm.GetCellMapping(cur_cell);
  const auto& adj_cell_mapping = m_sdm.GetCellMapping(adj_cell);

  const int i = cur_cell_mapping.MapFaceNode(ccf, ccfi);
  const auto& node_i_loc = cc_node_locs[i];

  const size_t adj_face_num_nodes = adj_cell_mapping.NumFaceNodes(acf);

  for (size_t fj=0; fj<adj_face_num_nodes; ++fj)
  {
    const int j = adj_cell_mapping.MapFaceNode(acf,fj);
    if ((node_i_loc - ac_node_locs[j]).NormSquare() < epsilon)
      return j;
  }

  throw std::logic_error(
    "lbs::acceleration::DiffusionMIPSolver::MapFaceNodeDisc: Mapping failure.");
}

double lbs::acceleration::DiffusionMIPSolver::
  CallLuaXYZFunction(lua_State* L, const std::string& lua_func_name,
                     const chi_mesh::Vector3& xyz)
{
  const std::string fname = "lbs::acceleration::DiffusionMIPSolver::"
                            "CallLuaXYZFunction";
  //============= Load lua function
  lua_getglobal(L, lua_func_name.c_str());

  //============= Error check lua function
  if (not lua_isfunction(L, -1))
    throw std::logic_error(fname + " attempted to access lua-function, " +
                           lua_func_name + ", but it seems the function"
                                           " could not be retrieved.");

  //============= Push arguments
  lua_pushnumber(L, xyz.x);
  lua_pushnumber(L, xyz.y);
  lua_pushnumber(L, xyz.z);

  //============= Call lua function
  //3 arguments, 1 result (double), 0=original error object
  double lua_return;
  if (lua_pcall(L,3,1,0) == 0)
  {
    LuaCheckNumberValue(fname, L, -1);
    lua_return = lua_tonumber(L,-1);
  }
  else
    throw std::logic_error(fname + " attempted to call lua-function, " +
                           lua_func_name + ", but the call failed." +
                           xyz.PrintStr());

  lua_pop(L,1); //pop the double, or error code

  return lua_return;
}