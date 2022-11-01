#ifndef FV_POLYGON_VALUES_H
#define FV_POLYGON_VALUES_H

#include "fv_cellbase.h"
#include "ChiMesh/Cell/cell.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//################################################################### Class def
namespace chi_math
{
  /**Finite Volume implementation for a polygon.
     *
     * - face_area[f] gives the area of the face.
     * - side_s_v[f][v] gives the vector for each leg of the
     *   triangle forming a side (face).*/
  class PolygonFVValues : public chi_math::CellFVValues
  {
  public:
    std::vector<std::vector<chi_mesh::Vector3>> side_legs;

    PolygonFVValues(const chi_mesh::Cell& poly_cell,
                    const chi_mesh::MeshContinuumConstPtr& grid);
  };
}

#endif