#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "chi_runtime.h"
#include "chi_log.h"

bool lbs::SteadySolver::
CheckPointInsideCell(const chi_mesh::Cell &cell,
                     const chi_mesh::MeshContinuum &grid_ref,
                     const chi_mesh::Vector3 &point)
{
  typedef chi_mesh::Vector3 Vec3;
  auto InsideTet = [](const Vec3& point,
                      const Vec3& v0, const Vec3& v1, const Vec3& v2)
  {
    const auto& v01 = v1-v0;
    const auto& v02 = v2-v0;

    const auto n = v01.Cross(v02).Normalized();
    const auto c = (v0 + v1 + v2)/3.0;

    const auto pc = point - c;

    if (pc.Dot(n) > 0.0)
      return true;
    else
      return false;
  };

  bool inside = true;
  if (cell.Type() == chi_mesh::CellType::SLAB)
  {
    const auto& v0 = grid_ref.vertices[cell.vertex_ids[0]];
    const auto& v1 = grid_ref.vertices[cell.vertex_ids[1]];

    const auto v01 = v1-v0;
    const auto v0p = point-v0;

    const double v0p_dot_v01 = v0p.Dot(v01);

    if (not (v0p_dot_v01 >= 0 and v0p_dot_v01 < v01.Norm()))
      inside = false;
  }//slab

  else if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    for (const auto& face : cell.faces)
    {
      const auto& vcp = point - face.centroid;

      if (vcp.Dot(face.normal) > 0)
      {
        inside = false;
        break;
      }
    }//for face
  }//polygon

  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    inside = false;
    //form tetra hedrons
    const auto& vcc = cell.centroid;
    for (const auto& face : cell.faces)
    {
      const auto& vfc = face.centroid;

      const size_t num_sides = face.vertex_ids.size();
      for (size_t s=0; s<num_sides; ++s)
      {
        const size_t sp1 = (s<(num_sides-1))? s+1 : 0;
        const auto& v0 = grid_ref.vertices[face.vertex_ids[s]];
        const auto& v1 = vfc;
        const auto& v2 = grid_ref.vertices[face.vertex_ids[sp1]];
        const auto& v3 = vcc;

        typedef std::tuple<Vec3, Vec3, Vec3> TetFace;

        std::vector<TetFace> tet_faces;
        tet_faces.emplace_back(v0, v1, v2);
        tet_faces.emplace_back(v0, v2, v3);
        tet_faces.emplace_back(v1, v3, v2);
        tet_faces.emplace_back(v0, v3, v1);

        bool inside_tet = true;
        for (const auto& tet_face : tet_faces)
        {
          if (not InsideTet(point, std::get<0>(tet_face),
                                   std::get<1>(tet_face),
                                   std::get<2>(tet_face)))
          {
            inside_tet = false;
            break;
          }
        }//for triangular tet_face
        if (inside_tet)
        {
          inside = true;
          break;
        }
      }//for side
      if (inside) break;
    }//for face
  }//polyhedron
  else
    return false; //TODO: Pop on no supported cell

  return inside;
}

void lbs::SteadySolver::InitializePointSources()
{
  const std::string fname = __FUNCTION__;

  typedef chi_math::SpatialDiscretization_PWLD PWLD;
  const auto& pwld = std::dynamic_pointer_cast<PWLD>(discretization);

  //============================================= Loop over point sources
  for (auto& point_source : point_sources)
  {
    if (point_source.Strength().size() != num_groups)
      throw std::logic_error(
        fname + ": Point source multigroup strength vector "
                "is not compatible with the number of "
                "groups in the simulaation. Expected " +
                std::to_string(num_groups) + " found " +
                std::to_string(point_source.Strength().size()));

    const auto& p = point_source.Location();
    double v_total = 0.0; //Total volume of all cells sharing
                          // this source
    std::vector<PointSource::ContainingCellInfo> temp_list;
    for (const auto& cell : grid->local_cells)
    {
      if (CheckPointInsideCell(cell, *grid, p))
      {
//        const auto& cell_view = pwld->GetCellMappingFE(cell.local_id);
        const auto& cell_view = pwld->GetCellMapping(cell);
        const auto& fe_values = pwld->GetUnitIntegrals(cell);
        const auto& M = fe_values.GetIntV_shapeI_shapeJ();
        const auto& I = fe_values.GetIntV_shapeI();

        std::vector<double> shape_values;
        cell_view.ShapeValues(point_source.Location(),
                              shape_values/**ByRef*/);

        const auto M_inv = chi_math::Inverse(M);

        const auto q_p_weights = chi_math::MatMul(M_inv, shape_values);

        double v_cell = 0.0;
        for (double val : I) v_cell += val;
        v_total += v_cell;

        temp_list.push_back(
          PointSource::ContainingCellInfo{v_cell,
                                          cell.local_id,
                                          shape_values,
                                          q_p_weights});
      }//if inside
    }//for local cell

    auto ghost_global_ids = grid->cells.GetGhostGlobalIDs();
    for (uint64_t ghost_global_id : ghost_global_ids)
    {
      const auto& neighbor_cell = grid->cells[ghost_global_id];
      if (CheckPointInsideCell(neighbor_cell, *grid, p))
      {
        const auto& neighbor_fe_values =
          pwld->GetUnitIntegrals(neighbor_cell);
        for (double val : neighbor_fe_values.GetIntV_shapeI())
          v_total += val;
      }//if point inside
    }//for ghost cell

    point_source.ClearInitializedInfo();
    for (const auto& info : temp_list)
    {
      point_source.AddContainingCellInfo(info.volume_weight/v_total,
                                         info.cell_local_id,
                                         info.shape_values,
                                         info.node_weights);
      const auto& cell = grid->local_cells[info.cell_local_id];
      //Output message
      {
        std::stringstream output;
        output << "Point source at " << p.PrintStr() << " assigned to cell "
               << cell.global_id << " with shape values ";
        for (double val : info.shape_values) output << val << " ";
        output << "volume_weight=" << info.volume_weight/v_total;

        chi::log.LogAll() << output.str();
      }
    }//for info in temp list
  }//for point_source
}