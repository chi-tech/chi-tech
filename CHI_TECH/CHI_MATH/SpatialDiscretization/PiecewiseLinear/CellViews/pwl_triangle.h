#ifndef _cfe_triangle_h
#define _cfe_triangle_h

#include "../pwl.h"
#include <vector>
#include "../../../Quadratures/quadrature.h"
#include "../../../../CHI_MESH/CHI_CELL/cell_triangle.h"


class TriangleFEView : public CellFEView
{
public:
  std::vector<chi_math::QuadraturePointXY*> qpoints;
  std::vector<double> w;
  double detJ;
  chi_mesh::Vector v01;
  chi_mesh::Vector v02;

  TriangleFEView(chi_mesh::CellTriangle* tri_cell,
                 chi_mesh::MeshContinuum* vol_continuum) : CellFEView(3)
  {
    //=========================================== Create the single quad point
    chi_math::QuadraturePointXY* new_qpoint = new chi_math::QuadraturePointXY;
    new_qpoint->x = 1.0/6.0;
    new_qpoint->y = 1.0/6.0;
    qpoints.push_back(new_qpoint);
    w.push_back(1.0/6.0);

    new_qpoint = new chi_math::QuadraturePointXY;
    new_qpoint->x = 4.0/6.0;
    new_qpoint->y = 1.0/6.0;
    qpoints.push_back(new_qpoint);
    w.push_back(1.0/6.0);

    new_qpoint = new chi_math::QuadraturePointXY;
    new_qpoint->x = 1.0/6.0;
    new_qpoint->y = 4.0/6.0;
    qpoints.push_back(new_qpoint);
    w.push_back(1.0/6.0);

    //=========================================== Get the vertices
    chi_mesh::Vertex v[3];
    for (int k=0;k<3;k++)
    {
      v[k] = *vol_continuum->nodes.at(tri_cell->v_index[k]);
    }

    //=========================================== Calculate legs
    v01 = v[1]-v[0];
    v02 = v[2]-v[0];

    //=========================================== Calculate determinant
    detJ = ((v01.x)*(v02.y) - (v02.x)*(v01.y));


  }

  double Varphi(int index, int qpoint_index)
  {
    if (index==0)
    {
      chi_math::QuadraturePointXY* qpoint = qpoints.at(qpoint_index);
      return 1.0 - qpoint->x - qpoint->y;
    }
    if (index==1)
    {
      chi_math::QuadraturePointXY* qpoint = qpoints.at(qpoint_index);
      return qpoint->x;
    }
    if (index==2)
    {
      chi_math::QuadraturePointXY* qpoint = qpoints.at(qpoint_index);
      return qpoint->y;
    }
    return 0;
  }

  double GradVarphi_x(int index, int qpoint_index)
  {

    if (index==0)
    {

      return (v02.y*(-1) - v01.y*(-1))/detJ;
    }
    if (index==1)
    {

      return (v02.y*(1) - v01.y*(0))/detJ;
    }
    if (index==2)
    {

      return (v02.y*(0) - v01.y*(1))/detJ;
    }

    return 0;
  }

  double GradVarphi_y(int index, int qpoint_index)
  {

    if (index==0)
    {

      return (-v02.x*(-1) + v01.x*(-1))/detJ;
    }
    if (index==1)
    {

      return (-v02.x*(1) + v01.x*(0))/detJ;
    }
    if (index==2)
    {

      return (-v02.x*(0) + v01.x*(1))/detJ;
    }


    return 0;
  }

  double DetJ(int qpoint_index)
  {
    return detJ;
  }
};

#endif