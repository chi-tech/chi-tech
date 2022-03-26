#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "chi_log.h"
extern ChiLog& chi_log;

void lbs::SteadySolver::InitializePointSources()
{
  const std::string fname = __FUNCTION__;
  typedef chi_mesh::Vector3 Vec3;
  auto InsideTet = [](const Vec3& point,
                      const Vec3& v0, const Vec3& v1, const Vec3& v2)
  {
    const auto& v01 = v1-v0;
    const auto& v02 = v2-v0;

    const auto n = v01.Cross(v02).Normalized();
    const auto c = (v0 + v1 + v2)/3.0;

    const auto pc = point - c;

    if (pc.Dot(n) < 0.0)
      return true;
    else
      return false;
  };

  typedef SpatialDiscretization_PWLD PWLD;
  const auto& pwld = std::dynamic_pointer_cast<PWLD>(discretization);

  for (auto& point_source : point_sources)
  {
    if (point_source.Strength().size() != num_groups)
      throw std::logic_error(fname + ": Point source multigroup strength vector "
                             "is not compatible with the number of "
                             "groups in the simulaation. Expected " +
                             std::to_string(num_groups) + " found " +
                             std::to_string(point_source.Strength().size()));

    const auto& p = point_source.Location();
    for (const auto& cell : grid->local_cells)
    {
      bool inside = true;
      if (cell.Type() == chi_mesh::CellType::SLAB)
      {
        const auto& v0 = grid->vertices[cell.vertex_ids[0]];
        const auto& v1 = grid->vertices[cell.vertex_ids[1]];

        const auto v01 = v1-v0;
        const auto v0p = p-v0;

        const double v0p_dot_v01 = v0p.Dot(v01);

        if (not (v0p_dot_v01 >= 0 and v0p_dot_v01 < v01.Norm()))
          inside = false;
      }//slab

      if (cell.Type() == chi_mesh::CellType::POLYGON)
      {
        for (const auto& face : cell.faces)
        {
          const auto& vcp = p - face.centroid;

          if (vcp.Dot(face.normal) > 0)
          {
            inside = false;
            break;
          }
        }//for face
      }//polygon

      if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
      {
        //form tetra hedrons
        const auto& vcc = cell.centroid;
        for (const auto& face : cell.faces)
        {
          const auto& vfc = face.centroid;

          const size_t num_sides = face.vertex_ids.size();
          for (size_t s=0; s<num_sides; ++s)
          {
            const size_t sp1 = (s<(num_sides-1))? s+1 : 0;
            const auto& v0 = grid->vertices[face.vertex_ids[s]];
            const auto& v1 = vfc;
            const auto& v2 = grid->vertices[face.vertex_ids[sp1]];
            const auto& v3 = vcc;

            typedef std::tuple<Vec3, Vec3, Vec3> TetFace;

            std::vector<TetFace> tet_faces;
            tet_faces.emplace_back(v0, v1, v2);
            tet_faces.emplace_back(v0, v2, v3);
            tet_faces.emplace_back(v1, v3, v2);
            tet_faces.emplace_back(v0, v3, v1);

            for (const auto& tet_face : tet_faces)
            {
              if (not InsideTet(p, std::get<0>(tet_face),
                                   std::get<1>(tet_face),
                                   std::get<2>(tet_face)))
              {
                inside = false;
                break;
              }
            }//for tri_face
            if (not inside) break;
          }//for side
          if (not inside) break;
        }//for face
      }//polyhedron

      if (inside)
      {
        const auto& cell_view = pwld->GetCellMappingFE(cell.local_id);

        std::vector<double> nodal_weights;
        cell_view->ShapeValues(point_source.Location(), nodal_weights/**ByRef*/);

        std::stringstream output;
        output << "Point source at " << p.PrintStr() << " assigned to cell "
               << cell.global_id << " with weights ";
        for (double val : nodal_weights)
          output << val << " ";

        chi_log.Log(LOG_ALL) << output.str();

        point_source.SetOwningCellLocalIDAndWeights(cell.local_id, nodal_weights);
        break;
      }
    }//for cell
  }//for point_source
}