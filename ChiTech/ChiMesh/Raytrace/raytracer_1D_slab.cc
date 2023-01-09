#include "raytracing.h"

#include "ChiMesh/Cell/cell.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Performs raytracing within a 1D-slab.*/
void chi_mesh::RayTracer::TraceSlab(const Cell &cell,
                                    Vector3 &pos_i,
                                    Vector3 &omega_i,
                                    bool& intersection_found,
                                    bool& backward_tolerance_hit,
                                    RayTracerOutputInformation& oi)
{
  const auto& grid = Grid();
  chi_mesh::Vector3 intersection_point;
  std::pair<double,double> weights;

  const double fabs_mu = std::fabs(omega_i.Dot(cell.faces[0].normal));

  double d_extend = (fabs_mu<1.0e-15)? 1.0e15 : extension_distance/fabs_mu;

  chi_mesh::Vector3 pos_f_line = pos_i + omega_i * d_extend;

  int num_faces = 2;
  for (int f=0; f<num_faces; f++)
  {
    uint64_t fpi = cell.vertex_ids[f]; //face point index
    chi_mesh::Vertex face_point = grid.vertices[fpi];

    bool intersects = chi_mesh::CheckPlaneLineIntersect(
      cell.faces[f].normal, face_point,
      pos_i, pos_f_line,
      intersection_point, &weights);

    double D = weights.first*d_extend;

    if ( (D > backward_tolerance) and intersects )
    {
      oi.distance_to_surface = D;
      oi.pos_f = intersection_point;

      oi.destination_face_index = f;
      oi.destination_face_neighbor = cell.faces[f].neighbor_id;
      intersection_found = true;
      break;
    }
    if (intersects)
      backward_tolerance_hit = true;
  }//for faces
}