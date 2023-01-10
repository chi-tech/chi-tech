#include "raytracing.h"

#include "ChiMesh/Cell/cell.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Performs raytracing within a 2D Polygon.*/
void chi_mesh::RayTracer::TracePolygon(const Cell &cell,
                                       Vector3 &pos_i,
                                       Vector3 &omega_i,
                                       bool& intersection_found,
                                       bool& backward_tolerance_hit,
                                       RayTracerOutputInformation& oi)
{
  const auto& grid = Grid();
  chi_mesh::Vector3 ip; //intersection point

  const double fabs_mu = std::fabs(omega_i.Dot(cell.faces[0].normal));

  double d_extend = (fabs_mu<1.0e-15)? 1.0e15 : extension_distance/fabs_mu;

  chi_mesh::Vector3 pos_f_line = pos_i + omega_i * d_extend;

  std::vector<RayTracerOutputInformation> face_intersections;

  size_t num_faces = cell.faces.size();
  face_intersections.reserve(num_faces);
  for (int f=0; f<num_faces; f++)
  {
    if (cell.faces[f].normal.Dot(omega_i) < 0.0) continue;

    RayTracerOutputInformation face_oi;

    uint64_t fpi = cell.faces[f].vertex_ids[0]; //face point index 0
    uint64_t fpf = cell.faces[f].vertex_ids[1]; //face point index 1
    const chi_mesh::Vertex& face_point_i = grid.vertices[fpi];
    const chi_mesh::Vertex& face_point_f = grid.vertices[fpf];

    bool intersects = chi_mesh::CheckLineIntersectStrip(
      face_point_i, face_point_f, cell.faces[f].normal,
      pos_i, pos_f_line, ip);

    double D = (ip - pos_i).Norm();

    if ( (D > backward_tolerance) and intersects )
    {
      face_oi.distance_to_surface = D;
      face_oi.pos_f = ip;

      face_oi.destination_face_index = f;
      face_oi.destination_face_neighbor = cell.faces[f].neighbor_id;
      intersection_found = true;
      face_intersections.emplace_back(std::move(face_oi));
      if (not perform_concavity_checks) break;
    }//if intersects
    if ( (D < backward_tolerance) and intersects )
      backward_tolerance_hit = true;
  }//for faces

  //======================================== Determine closest intersection
  if (not perform_concavity_checks and not face_intersections.empty())
    oi = face_intersections.back();
  else if (perform_concavity_checks and not face_intersections.empty())
  {
    auto closest_intersection = &face_intersections.back();
    for (auto& intersection : face_intersections)
      if (intersection.distance_to_surface <
          closest_intersection->distance_to_surface)
        closest_intersection = &intersection;

    oi = *closest_intersection;
  }
  else
  {
    RayTracerOutputInformation blank_oi;
    blank_oi.distance_to_surface = 1.0e15;
    oi = blank_oi;
  }
}