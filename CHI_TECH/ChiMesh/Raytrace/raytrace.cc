#include "raytracing.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>

#include <chi_log.h>

extern ChiLog& chi_log;

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
  const chi_mesh::Vector3 &pos_i,
  const chi_mesh::Vector3 &omega_i,
  double& d_to_surface,
  chi_mesh::Vector3 &pos_f,
  double epsilon_nudge,
  double backward_tolerance,
  double extension_distance,
  int func_depth)
{
  chi_mesh::RayDestinationInfo dest_info;

//  double extention_distance = 1.0e5;
  extension_distance = d_to_surface;
  chi_mesh::Vector3 pos_f_line = pos_i + omega_i * extension_distance;

  bool intersection_found = false;
  bool backward_tolerance_hit = false;

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
  if (cell->Type() == chi_mesh::CellType::SLAB)
  {
    auto slab_cell = (chi_mesh::CellSlab*)cell;

    chi_mesh::Vector3 intersection_point;
    std::pair<double,double> weights;

    int num_faces = 2;
    for (int f=0; f<num_faces; f++)
    {
      int fpi = slab_cell->vertex_ids[f]; //face point index
      chi_mesh::Vertex face_point = *grid->vertices[fpi];

      bool intersects = chi_mesh::CheckPlaneLineIntersect(
        slab_cell->faces[f].normal, face_point,
        pos_i, pos_f_line,
        intersection_point, weights);

      double D = weights.first*extension_distance;

      if ( (D > backward_tolerance) and intersects )
      {
        d_to_surface = D;
        pos_f = intersection_point;

        dest_info.destination_face_index = f;
        dest_info.destination_face_neighbor = slab_cell->faces[f].neighbor;
        intersection_found = true;
        break;
      }
      if (intersects)
        backward_tolerance_hit = true;
    }//for faces
  }//slab
  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ POLYGON
  else if (cell->Type() == chi_mesh::CellType::POLYGON)
  {
    auto poly_cell = (chi_mesh::CellPolygon*)cell;

    chi_mesh::Vector3 ip; //intersetion point

    int num_faces = poly_cell->faces.size();
    for (int f=0; f<num_faces; f++)
    {
      int fpi = poly_cell->faces[f].vertex_ids[0]; //face point index 0
      int fpf = poly_cell->faces[f].vertex_ids[1]; //face point index 1
      chi_mesh::Vertex& face_point_i = *grid->vertices[fpi];
      chi_mesh::Vertex& face_point_f = *grid->vertices[fpf];

      bool intersects = chi_mesh::CheckLineIntersectStrip(
        face_point_i, face_point_f, poly_cell->faces[f].normal,
        pos_i, pos_f_line, ip);

      double D = (ip - pos_i).Norm();

      if ( (D > backward_tolerance) and intersects )
      {
        d_to_surface = D;
        pos_f = ip;

        dest_info.destination_face_index = f;
        dest_info.destination_face_neighbor = poly_cell->faces[f].neighbor;
        intersection_found = true;
        break;
      }//if intersects
      if (intersects)
        backward_tolerance_hit = true;
    }//for faces
  }//polygon
  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ POLYHEDRON
  else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
  {
    auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;

    chi_mesh::Vector3 ip = pos_i; //Intersection point

    int num_faces = polyh_cell->faces.size();
    for (int f=0; f<num_faces; f++)
    {
      auto& edges = polyh_cell->GetFaceEdges(f);

      size_t num_sides = edges.size();
      for (size_t s=0; s<num_sides; s++)
      {
        chi_mesh::Vertex& v0 = *grid->vertices[edges[s][0]];
        chi_mesh::Vertex& v1 = *grid->vertices[edges[s][1]];
        chi_mesh::Vertex& v2 = polyh_cell->faces[f].centroid;

        bool intersects =
          chi_mesh::CheckLineIntersectTriangle2(v0,v1,v2,pos_i,omega_i,ip);

        if (intersects)
        {
          d_to_surface = (ip - pos_i).Norm();
          pos_f = ip;

          dest_info.destination_face_index = f;
          dest_info.destination_face_neighbor = polyh_cell->faces[f].neighbor;

          intersection_found = true;
          break;
        }//if intersects
      }//for side

      if (intersection_found) break;
    }//for faces
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
    if (func_depth < 5)
    {
      //chi_log.Log(LOG_ALLERROR) << "Particle nudged";
      //Vector from position to cell-centroid
      chi_mesh::Vector3 v_p_i_cc = (cell->centroid - pos_i);
      chi_mesh::Vector3 pos_i_nudged = pos_i + v_p_i_cc * epsilon_nudge;

      //printf("%.12f %.12f %.12f\n",pos_i_nudged.x,pos_i_nudged.y,pos_i_nudged.z);

      dest_info =
        RayTrace(grid,cell,pos_i_nudged,omega_i,d_to_surface,pos_f,
                 epsilon_nudge,backward_tolerance,extension_distance,
                 func_depth+1);

      return dest_info;
    }

    if (func_depth < 7)
    {
      //chi_log.Log(LOG_ALLERROR) << "Particle nudged";
      //Vector from position to cell-centroid
      chi_mesh::Vector3 v_p_i_cc = (cell->centroid - pos_i).Cross(omega_i);
      chi_mesh::Vector3 pos_i_nudged = pos_i + v_p_i_cc * epsilon_nudge;

      //printf("%.12f %.12f %.12f\n",pos_i_nudged.x,pos_i_nudged.y,pos_i_nudged.z);

      dest_info =
        RayTrace(grid,cell,pos_i_nudged,omega_i,d_to_surface,pos_f,
                 epsilon_nudge,backward_tolerance,extension_distance,
                 func_depth+1);

      return dest_info;
    }


    std::stringstream outstr;

    outstr
      << "Intersection not found at function level " << func_depth << "."
      << ((backward_tolerance_hit)? " Backward tolerance hit. " : "")
      << "For particle xyz="
      << pos_i.PrintS() << " uvw="
      << omega_i.PrintS() << " in cell " << cell->global_id
      << " with vertices: \n";

    for (auto vi : cell->vertex_ids)
      outstr << grid->vertices[vi]->PrintS() << "\n";

    int f=0;
    for (auto& face : cell->faces)
    {
      outstr << "Face with centroid: " << face.centroid.PrintS() << " ";
      outstr << "n=" << face.normal.PrintS() << "\n";
      for (auto vi : face.vertex_ids)
        outstr << grid->vertices[vi]->PrintS() << "\n";

      //TODO: Temp code
      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;

      auto& edges = polyh_cell->GetFaceEdges(f);
      for (auto& edge : edges)
        outstr << "Edge " << grid->vertices[edge[0]]->PrintS() << " "
                          << grid->vertices[edge[1]]->PrintS() << "\n";

               ++f;
      //TODO: Temp code
    }



    chi_log.Log(LOG_ALLERROR) << outstr.str();
    exit(EXIT_FAILURE);
  }

  return dest_info;
}
