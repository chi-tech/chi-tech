#include "raytracing.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include <ChiMesh/Cell/cell_slabv2.h>
#include <ChiMesh/Cell/cell_polygonv2.h>
#include <ChiMesh/Cell/cell_polyhedronv2.h>

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**Performs raytracing on chi_mesh type cells. This algorithm essentially
 * looks for intersections with planes, strips or triangles.
 *
 * \param grid          Reference to the grid on which the specific cell resides.
 * \param cell          A pointer to the cell base.
 * \param pos_i         The initial position of the ray.
 * \param omega_i       The direction of the ray.
 * \param d_to_surface  A reference to a variable that will receive the distance to the
 *        surface.
 * \param pos_f         A reference to the vector to receive the final position.
 * \param get_segments  (Optional) A flag indicating .*/
chi_mesh::RayDestinationInfo chi_mesh::RayTrace(
  chi_mesh::MeshContinuum* grid,
  chi_mesh::Cell *cell,
  const chi_mesh::Vector &pos_i,
  const chi_mesh::Vector &omega_i,
  double& d_to_surface,
  chi_mesh::Vector &pos_f)
{
  chi_mesh::RayDestinationInfo dest_info;

  double extention_distance = 1.0e15;
  chi_mesh::Vector pos_f_line = pos_i + omega_i*extention_distance;

  bool intersection_found = false;

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
  if (cell->Type() == chi_mesh::CellType::SLABV2)
  {
    auto slab_cell = (chi_mesh::CellSlabV2*)cell;

    chi_mesh::Vector intersection_point;
    std::pair<double,double> weights;

    int num_faces = 2;
    for (int f=0; f<num_faces; f++)
    {
      int fpi = slab_cell->vertex_ids[f]; //face point index
      chi_mesh::Vertex face_point = *grid->nodes[fpi];

      bool intersects = chi_mesh::CheckPlaneLineIntersect(
        slab_cell->faces[f].normal, face_point,
        pos_i, pos_f_line,
        intersection_point, weights);

      double D = weights.first*extention_distance;

      if ( (D > 1.0e-10) and intersects )
      {
        d_to_surface = D;
        pos_f = intersection_point;

        dest_info.destination_face_index = f;
        dest_info.destination_face_neighbor = slab_cell->faces[f].neighbor;
        intersection_found = true;
        break;
      }
    }//for faces


  }//slab
  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ POLYGON
  else if (cell->Type() == chi_mesh::CellType::POLYGONV2)
  {
    auto poly_cell = (chi_mesh::CellPolygonV2*)cell;

    chi_mesh::Vector intersection_point;

    int num_faces = poly_cell->faces.size();
    for (int f=0; f<num_faces; f++)
    {
      int fpi = poly_cell->faces[f].vertex_ids[0]; //face point index 0
      int fpf = poly_cell->faces[f].vertex_ids[1]; //face point index 1
      chi_mesh::Vertex face_point_i = *grid->nodes[fpi];
      chi_mesh::Vertex face_point_f = *grid->nodes[fpf];

      bool intersects = chi_mesh::CheckLineIntersectStrip(
        face_point_i,
        face_point_f,
        poly_cell->faces[f].normal,
        pos_i, pos_f_line,
        intersection_point);

      if (intersects)
      {
        double D = (intersection_point - pos_i).Norm();

        if (D > 1.0e-10)
        {
          d_to_surface = D;
          pos_f = intersection_point;

          dest_info.destination_face_index = f;
          dest_info.destination_face_neighbor = poly_cell->faces[f].neighbor;
          intersection_found = true;
          break;
        }
      }//if intersects
    }//for faces
  }//polygon
  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ POLYHEDRON
  else if (cell->Type() == chi_mesh::CellType::POLYHEDRONV2)
  {
    auto polyh_cell = (chi_mesh::CellPolyhedronV2*)cell;

    chi_mesh::Vector intersection_point;

    int num_faces = polyh_cell->faces.size();
    for (int f=0; f<num_faces; f++)
    {
      auto& edges = polyh_cell->GetFaceEdges(f);

      size_t num_sides = edges.size();
      for (size_t s=0; s<num_sides; s++)
      {
        int v0i = edges[s][0]; //face point index 0
        int v1i = edges[s][1]; //face point index 1
        chi_mesh::Vertex& v0 = *grid->nodes[v0i];
        chi_mesh::Vertex& v1 = *grid->nodes[v1i];
        chi_mesh::Vertex& v2 = polyh_cell->faces[f].centroid;

        auto v01 = v1 - v0;
        auto v12 = v2 - v1;

        auto n = v01.Cross(v12);
        n = n/n.Norm();

        std::pair<double,double> weights;

        bool intersects = chi_mesh::CheckPlaneLineIntersect(
          n, v0, pos_i, pos_f_line,
          intersection_point, weights);



        if (intersects)
        {
//          chi_log.Log(LOG_0) << intersection_point.PrintS();
//          chi_log.Log(LOG_0) << "    v0: " << v0.PrintS();
//          chi_log.Log(LOG_0) << "    v1: " << v1.PrintS();
//          chi_log.Log(LOG_0) << "    v2: " << v2.PrintS();

          intersects = chi_mesh::CheckPointInTriangle(
            v0,v1,v2,n,intersection_point);
        }


        if (intersects)
        {
          double D = (intersection_point - pos_i).Norm();

          if (D > 1.0e-10)
          {
            d_to_surface = D;
            pos_f = intersection_point;

            dest_info.destination_face_index = f;
            dest_info.destination_face_neighbor = polyh_cell->faces[f].neighbor;

//            chi_log.Log(LOG_0)
//              << "Distance to surface"
//            usleep(100000);
            intersection_found = true;
            break;
          }
        }//if intersects
      }


    }//for faces
    if (!intersection_found)
    {
      chi_log.Log(LOG_0) << "Cell: ";
      for (auto vi : polyh_cell->vertex_ids)
        chi_log.Log(LOG_0) << grid->nodes[vi]->PrintS();
    }
  }//polyhedron
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "Unsupported cell type encountered in call to "
      << "chi_mesh::RayTrace.";
    exit(EXIT_FAILURE);
  }

  if (!intersection_found)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Intersection not found. "
      << pos_i.PrintS() << " "
      << omega_i.PrintS();
    exit(EXIT_FAILURE);
    usleep(1000000);
  }

  return dest_info;
}
