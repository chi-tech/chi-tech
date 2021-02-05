#ifndef _fv_polyhedron_h
#define _fv_polyhedron_h

#include "fv_cellbase.h"
#include <ChiMesh/Cell/cell_polyhedron.h>
#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>

//################################################################### Class def
/**Finite Volume implementation for a polyhedron.
 *
 * - face_area[f] gives the area of the face.
 * - face_side_v[f][s][v] gives the vector for each leg of the
 *   tetrahedron forming a sides.*/
class PolyhedronFVView : public CellFVView
{
private:
  chi_mesh::MeshContinuumPtr grid;

public:
  std::vector<std::vector<double>>           face_side_area;
  std::vector<std::vector<double>>           face_side_volume;
  std::vector<std::vector<std::vector<chi_mesh::Vector3>>> face_side_vectors;

  PolyhedronFVView(chi_mesh::CellPolyhedron* polyh_cell,
                chi_mesh::MeshContinuumPtr vol_continuum) :
                CellFVView(polyh_cell->vertex_ids.size())
  {
    grid = vol_continuum;

    volume = 0.0;
    chi_mesh::Vector3& vcc = polyh_cell->centroid;

    int num_faces = polyh_cell->faces.size();
    face_side_vectors.resize(num_faces);
    face_side_area.resize(num_faces);
    face_side_volume.resize(num_faces);
    face_area.resize(num_faces,0.0);
    for (int f=0; f<num_faces; f++)
    {
      std::vector<std::vector<int>> edges = polyh_cell->GetFaceEdges(f);

      for (auto edge : edges)
      {
        int v0i = edge[0];
        int v1i = edge[1];

        chi_mesh::Vector3& v0 = *grid->vertices[v0i];
        chi_mesh::Vector3& v1 = polyh_cell->faces[f].centroid;
        chi_mesh::Vector3& v2 = *grid->vertices[v1i];
        chi_mesh::Vector3& v3 = vcc;

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