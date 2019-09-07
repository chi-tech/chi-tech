#include "pwl_polyhedron.h"

double PolyhedronFEView::TetShape(int index,
                                  int qpoint_index,
                                  bool on_surface)
{
  QPointXYZ* qpoint;
  if (!on_surface)
    qpoint = quadratures[DEG3]->qpoints.at(qpoint_index);
  else
    qpoint = quadratures[DEG3_SURFACE]->qpoints.at(qpoint_index);
  double value = 0.0;

  if (index == 0)
  {value = 1.0 - qpoint->x - qpoint->y - qpoint->z;}
  if (index == 1){value = qpoint->x;}
  if (index == 2){value = qpoint->y;}
  if (index == 3){value = qpoint->z;}

  return value;
}

double PolyhedronFEView::TetGradShape_x(int index, int qpoint_index)
{
  double value = 0.0;
  if (index == 0){value = -1.0;}
  if (index == 1){value =  1.0;}
  if (index == 2){value =  0.0;}
  if (index == 3){value =  0.0;}

  return value;
}

double PolyhedronFEView::TetGradShape_y(int index, int qpoint_index)
{
  double value = 0.0;
  if (index == 0){value = -1.0;}
  if (index == 1){value =  0.0;}
  if (index == 2){value =  1.0;}
  if (index == 3){value =  0.0;}

  return value;
}

double PolyhedronFEView::TetGradShape_z(int index, int qpoint_index)
{
  double value = 0.0;
  if (index == 0){value = -1.0;}
  if (index == 1){value =  0.0;}
  if (index == 2){value =  0.0;}
  if (index == 3){value =  1.0;}

  return value;
}