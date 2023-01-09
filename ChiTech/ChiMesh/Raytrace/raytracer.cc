#include "raytracing.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

const chi_mesh::MeshContinuum& chi_mesh::RayTracer::Grid() const
{
  return reference_grid;
}

//###################################################################
chi_mesh::RayTracerOutputInformation chi_mesh::RayTracer::
  TraceRay(const Cell &cell,
           Vector3 &pos_i,
           Vector3 &omega_i,
           int function_depth/*=0*/)
{
  if (not cell_sizes.empty())
    SetTolerancesFromCellSize(cell_sizes[cell.local_id]);

  RayTracerOutputInformation oi;

  bool intersection_found = false;
  bool backward_tolerance_hit = false;

  if (cell.Type() == chi_mesh::CellType::SLAB)
    TraceSlab(cell, pos_i, omega_i,
              intersection_found    /*byRef*/,
              backward_tolerance_hit/*byRef*/,
              oi                    /*byRef*/);
  else if (cell.Type() == chi_mesh::CellType::POLYGON)
    TracePolygon(cell, pos_i, omega_i,
                 intersection_found    /*byRef*/,
                 backward_tolerance_hit/*byRef*/,
                 oi                    /*byRef*/);
  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    TracePolyhedron(cell, pos_i, omega_i,
                    intersection_found    /*byRef*/,
                    backward_tolerance_hit/*byRef*/,
                    oi                    /*byRef*/);
  else
    throw std::logic_error("Unsupported cell type encountered in call to "
                           "chi_mesh::RayTrace.");

  if (!intersection_found)
  {
    if (function_depth < 5)
    {
      // Nudge particle towards centroid
      chi_mesh::Vector3 v_p_i_cc = (cell.centroid - pos_i);
      chi_mesh::Vector3 pos_i_nudged = pos_i + v_p_i_cc * epsilon_nudge;

      oi = TraceRay(cell,pos_i_nudged,omega_i,function_depth+1);

      return oi;
    }

    if (function_depth < 7)
    {
      // Nudge particle away from line between location and cell center
      chi_mesh::Vector3 v_p_i_cc = (cell.centroid - pos_i).Cross(omega_i);
      chi_mesh::Vector3 pos_i_nudged = pos_i + v_p_i_cc * epsilon_nudge;

      oi = TraceRay(cell,pos_i_nudged,omega_i,function_depth+1);

      return oi;
    }


    std::stringstream outstr;

    outstr
      << "Intersection not found at function level " << function_depth << "."
      << ((backward_tolerance_hit)? " Backward tolerance hit. " : "")
      << "For particle xyz="
      << pos_i.PrintS() << " uvw="
      << omega_i.PrintS() << " " << (pos_i + extension_distance*omega_i).PrintS()
      << " " << extension_distance
      << " in cell " << cell.global_id
      << " with vertices: \n";

    const auto& grid = Grid();

    for (auto vi : cell.vertex_ids)
      outstr << grid.vertices[vi].PrintS() << "\n";

    for (auto& face : cell.faces)
    {
      outstr << "Face with centroid: " << face.centroid.PrintS() << " ";
      outstr << "n=" << face.normal.PrintS() << "\n";
      for (auto vi : face.vertex_ids)
        outstr << grid.vertices[vi].PrintS() << "\n";
    }

    outstr << "o Cell\n";
    for (auto& vid : cell.vertex_ids)
    {
      auto& v = grid.vertices[vid];
      outstr << "v " << v.x << " " << v.y << " " << v.z << "\n";
    }

    for (auto& face : cell.faces)
    {
      auto& v = face.centroid;
      outstr << "v " << v.x << " " << v.y << " " << v.z << "\n";
    }

    for (size_t f=0; f < cell.faces.size(); ++f)
    {
      auto& face = cell.faces[f];
      outstr << "f ";
      for (auto vid : face.vertex_ids)
      {
        size_t ref_cell_id=0;
        for (size_t cid=0; cid < cell.vertex_ids.size(); ++cid)
          if (cell.vertex_ids[cid] == vid)
            ref_cell_id = cid + 1;

        outstr << ref_cell_id << "// ";
      }
      outstr << "\n";
    }

    oi.particle_lost = true;
    oi.lost_particle_info = outstr.str();
  }

  return oi;
}


//###################################################################
chi_mesh::RayTracerOutputInformation chi_mesh::RayTracer::
TraceIncidentRay(const Cell& cell,
                 const Vector3& pos_i,
                 const Vector3& omega_i)
{
  const auto cell_type = cell.Type();
  const double cell_char_length = cell_sizes[cell.local_id];
  const auto& grid = reference_grid;

  bool intersects_cell = false;
  chi_mesh::Vector3 I;

  size_t f = 0;
  for (const auto& face : cell.faces)
  {
    if (face.normal.Dot(omega_i) > 0.0) { ++f; continue/*the loop*/; }

    const auto& p0 = grid.vertices[face.vertex_ids[0]];
    const auto& n = face.normal;

    const auto ppos_i = p0 - pos_i;
    const double d = ppos_i.Dot(omega_i);

    const auto pos_ext = pos_i + (d + cell_char_length)*omega_i;

    using namespace chi_mesh;
    {
      if (cell_type == CellType::SLAB)
      {
        intersects_cell = CheckPlaneLineIntersect(n, p0, pos_i, pos_ext, I);
      }//SLAB
      else if (cell_type == CellType::POLYGON)
      {
        const auto& p1 = grid.vertices[face.vertex_ids[1]];
        intersects_cell = CheckLineIntersectStrip(p0, p1, n, pos_i, pos_ext, I);
      }//POLYGON
      else if (cell_type == CellType::POLYHEDRON)
      {
        const auto& vids = face.vertex_ids;
        const size_t num_sides = face.vertex_ids.size();
        for (size_t s=0; s<num_sides; ++s)
        {
          uint64_t v0i = vids[s];
          uint64_t v1i = (s < (num_sides-1))? vids[s+1] : vids[0];

          const auto& v0 = grid.vertices[v0i];
          const auto& v1 = grid.vertices[v1i];
          const auto& v2 = face.centroid;


          const auto v01 = v1 - v0;
          const auto v02 = v2 - v0;
          const auto n_est = v01.Cross(v02);

          if (n_est.Dot(omega_i) > 0.0) continue;

          intersects_cell = CheckLineIntersectTriangle2(v0, v1, v2,
                                                        pos_i, omega_i, I);
          if (intersects_cell) break;
        }//for side
      }//POLYHEDRON
    }
    if (intersects_cell) break;
    ++f;
  }//for face

  RayTracerOutputInformation oi;
  if (intersects_cell)
  {
    oi.distance_to_surface       = (I - pos_i).Norm();
    oi.pos_f                     = I;
    oi.destination_face_index    = f;
    oi.destination_face_neighbor = cell.global_id;
    oi.particle_lost             = false;
  }
  else
    oi.particle_lost = true;

  return oi;
}
