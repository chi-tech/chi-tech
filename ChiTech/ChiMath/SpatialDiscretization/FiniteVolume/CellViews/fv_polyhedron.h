#ifndef FV_POLYHEDRON_VALUES_H
#define FV_POLYHEDRON_VALUES_H

#include "fv_cellbase.h"
#include "ChiMesh/Cell/cell.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//################################################################### Class def
namespace chi_math
{
  /**Finite Volume implementation for a polyhedron.
     *
     * - face_area[f] gives the area of the face.
     * - face_side_v[f][s][v] gives the vector for each leg of the
     *   tetrahedron forming a sides.*/
  class PolyhedronFVValues : public chi_math::CellFVValues
  {
  public:
    std::vector<std::vector<double>>           face_side_area;
    std::vector<std::vector<double>>           face_side_volume;
    std::vector<std::vector<std::vector<chi_mesh::Vector3>>> face_side_vectors;

    PolyhedronFVValues(const chi_mesh::Cell& polyh_cell,
                       const chi_mesh::MeshContinuumConstPtr& grid);
  };
}

#endif