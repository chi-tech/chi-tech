#include "pwl_polyhedron.h"

/**Precomputes the shape function values of a face-side pair
 * at a quadrature point*/
double PolyhedronPWLFEValues::FaceSideShape(int face_index, int side_index,
                                            int i, int qpoint_index, bool on_surface)
{
  double value = 0.0;
  int    index = node_side_maps[i].face_map[face_index].
    side_map[side_index].index;
  double betaf = face_betaf[face_index];

  value += TetShape(index, qpoint_index, on_surface);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    value += betaf* TetShape(1, qpoint_index, on_surface);}
  value += alphac* TetShape(3, qpoint_index, on_surface);

  return value;
}

/**Precomputes the gradx-shape function values of a face-side pair
 * at a quadrature point*/
double PolyhedronPWLFEValues::FaceSideGradShape_x(int face_index,
                                                  int side_index,
                                                  int i)
{
  double value = 0.0;
  double tetdfdx = 0.0;
  double tetdfdy = 0.0;
  double tetdfdz = 0.0;
  int    index = node_side_maps[i].face_map[face_index].
    side_map[side_index].index;
  double betaf = face_betaf[face_index];

  tetdfdx += TetGradShape_x(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdx += betaf* TetGradShape_x(1);}
  tetdfdx += alphac* TetGradShape_x(3);

  tetdfdy += TetGradShape_y(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdy += betaf* TetGradShape_y(1);}
  tetdfdy += alphac* TetGradShape_y(3);

  tetdfdz += TetGradShape_z(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdz += betaf* TetGradShape_z(1);}
  tetdfdz += alphac* TetGradShape_z(3);

  value += face_data[face_index].sides[side_index].JTinv.GetIJ(0, 0) * tetdfdx;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(0, 1) * tetdfdy;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(0, 2) * tetdfdz;

  return value;
}

/**Precomputes the grady-shape function values of a face-side pair
 * at a quadrature point*/
double PolyhedronPWLFEValues::FaceSideGradShape_y(int face_index,
                                                  int side_index,
                                                  int i)
{
  double value = 0.0;
  double tetdfdx = 0.0;
  double tetdfdy = 0.0;
  double tetdfdz = 0.0;
  int    index = node_side_maps[i].face_map[face_index].
    side_map[side_index].index;
  double betaf = face_betaf[face_index];

  tetdfdx += TetGradShape_x(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdx += betaf* TetGradShape_x(1);}
  tetdfdx += alphac* TetGradShape_x(3);

  tetdfdy += TetGradShape_y(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdy += betaf* TetGradShape_y(1);}
  tetdfdy += alphac* TetGradShape_y(3);

  tetdfdz += TetGradShape_z(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdz += betaf* TetGradShape_z(1);}
  tetdfdz += alphac* TetGradShape_z(3);

  value += face_data[face_index].sides[side_index].JTinv.GetIJ(1, 0) * tetdfdx;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(1, 1) * tetdfdy;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(1, 2) * tetdfdz;

  return value;
}

/**Precomputes the gradz-shape function values of a face-side pair
 * at a quadrature point*/
double PolyhedronPWLFEValues::FaceSideGradShape_z(int face_index,
                                                  int side_index,
                                                  int i)
{
  double value = 0.0;
  double tetdfdx = 0.0;
  double tetdfdy = 0.0;
  double tetdfdz = 0.0;
  int    index = node_side_maps[i].face_map[face_index].
    side_map[side_index].index;
  double betaf = face_betaf[face_index];

  tetdfdx += TetGradShape_x(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdx += betaf* TetGradShape_x(1);}
  tetdfdx += alphac* TetGradShape_x(3);

  tetdfdy += TetGradShape_y(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdy += betaf* TetGradShape_y(1);}
  tetdfdy += alphac* TetGradShape_y(3);

  tetdfdz += TetGradShape_z(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdz += betaf* TetGradShape_z(1);}
  tetdfdz += alphac* TetGradShape_z(3);

  value += face_data[face_index].sides[side_index].JTinv.GetIJ(2, 0) * tetdfdx;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(2, 1) * tetdfdy;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(2, 2) * tetdfdz;

  return value;
}