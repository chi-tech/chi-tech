#include "fv_polygon.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

chi_math::PolygonFVValues::
  PolygonFVValues(const chi_mesh::Cell& poly_cell,
                  const chi_mesh::MeshContinuumConstPtr& grid) :
  CellFVValues(grid,
               poly_cell.centroid,
               std::vector<std::vector<int>>(poly_cell.faces.size(),{-1}))
{
  volume = 0.0;

  size_t num_faces = poly_cell.faces.size();
  side_legs.resize(num_faces);
  for (size_t f=0; f<num_faces; f++)
  {
    uint64_t v0i = poly_cell.faces[f].vertex_ids[0];
    uint64_t v1i = poly_cell.faces[f].vertex_ids[1];

    const auto& v0 = grid->vertices[v0i];
    const auto& v1 = grid->vertices[v1i];
    chi_mesh::Vector3 v2 = poly_cell.centroid;

    face_area.push_back((v1-v0).Norm());

    chi_mesh::Vector3 sidev01 = v1 - v0;
    chi_mesh::Vector3 sidev02 = v2 - v0;

    side_legs[f].push_back(sidev01);
    side_legs[f].push_back(sidev02);

    double sidedetJ = ((sidev01.x)*(sidev02.y) - (sidev02.x)*(sidev01.y));

    volume += sidedetJ/2.0;
  }

}