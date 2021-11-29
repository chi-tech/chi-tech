#ifndef FV_POLYHEDRON_VALUES_H
#define FV_POLYHEDRON_VALUES_H

#include "fv_cellbase.h"
#include "ChiMesh/Cell/cell.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//################################################################### Class def
/**Finite Volume implementation for a polyhedron.
 *
 * - face_area[f] gives the area of the face.
 * - face_side_v[f][s][v] gives the vector for each leg of the
 *   tetrahedron forming a sides.*/
class PolyhedronFVValues : public CellFVValues
{
public:
  std::vector<std::vector<double>>           face_side_area;
  std::vector<std::vector<double>>           face_side_volume;
  std::vector<std::vector<std::vector<chi_mesh::Vector3>>> face_side_vectors;

  PolyhedronFVValues(const chi_mesh::Cell& polyh_cell,
                     const chi_mesh::MeshContinuum& grid) :
    CellFVValues(polyh_cell.vertex_ids.size())
  {
    volume = 0.0;
    auto& vcc = polyh_cell.centroid;

    size_t num_faces = polyh_cell.faces.size();
    face_side_vectors.resize(num_faces);
    face_side_area.resize(num_faces);
    face_side_volume.resize(num_faces);
    face_area.resize(num_faces,0.0);
    for (size_t f=0; f<num_faces; f++)
    {
      const auto& face = polyh_cell.faces[f];
      const size_t num_edges = face.vertex_ids.size();
      for (size_t e=0; e<num_edges; ++e)
      {
        size_t ep1 = (e < (num_edges-1))? e+1 : 0;
        uint64_t v0i = face.vertex_ids[e  ];
        uint64_t v1i = face.vertex_ids[ep1];

        const auto& v0 = grid.vertices[v0i];
        const auto& v1 = polyh_cell.faces[f].centroid;
        const auto& v2 = grid.vertices[v1i];
        const auto& v3 = vcc;

        std::vector<chi_mesh::Vector3> side_legs(3);
        side_legs[0] = v1-v0;
        side_legs[1] = v2-v0;
        side_legs[2] = v3-v0;

        face_side_vectors[f].push_back(side_legs);

        chi_mesh::Matrix3x3 J;

        J.SetColJVec(0,side_legs[0]);
        J.SetColJVec(1,side_legs[1]);
        J.SetColJVec(2,side_legs[2]);

        chi_mesh::Vector3& sidev01 = side_legs[0];
        chi_mesh::Vector3& sidev02 = side_legs[1];

        double side_area = (sidev01.Cross(sidev02)).Norm()/2.0;
        face_side_area[f].emplace_back(side_area);
        face_area[f] += side_area;

        double side_volume = J.Det()/6.0;

        face_side_volume[f].emplace_back(side_volume);
        volume += side_volume;
      }//for edge

    }//for face
//    std::cout << volume << "\n";
  }
};

#endif