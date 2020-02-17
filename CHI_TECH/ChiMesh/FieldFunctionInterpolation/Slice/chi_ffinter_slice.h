#ifndef _chi_ffinter_slice_h
#define _chi_ffinter_slice_h

#include "../chi_ffinterpolation.h"
#include "../../chi_mesh.h"

#include <petscksp.h>

#define FFI_PROP_SLICEPOINT 1
#define FFI_PROP_SLICENORMAL 2
#define FFI_PROP_SLICETANGENT 3
#define FFI_PROP_SLICEBINORM 4
#define FFI_PROP_OPERATION 5
#define FFI_PROP_LOGICAL_VOLUME 8

struct FFIFaceEdgeIntersection
{
  int v0_g_index, v1_g_index;
  int v0_dofindex_cell, v1_dofindex_cell;
  chi_mesh::Vector3 point;
  chi_mesh::Vector3 point2d;
  std::pair<double,double> weights;
  double point_value;

  FFIFaceEdgeIntersection()
  {
    point_value = 0.0;
  }
};

struct FFICellIntersection
{
  int cell_global_index;
  std::vector<FFIFaceEdgeIntersection*> intersections;
  chi_mesh::Vector3 intersection_centre;
  chi_mesh::Vector3 intersection_2d_centre;
  double cell_avg_value;

  FFICellIntersection()
  {
    cell_avg_value = 0.0;
  }
};

//###################################################################
/** A slice based interpolation function.
 *
 * This functionality needs to cater for numerous spatial discretizations.
 * The most simple one is cell-averaged values and the more complicated ones
 * are PWLD and then CFEM.
 *
 * Cell average values requires computing the slice of the polyhedron and then
 * computing the centroid of that cut. This can be done cell by cell.*/
class chi_mesh::FieldFunctionInterpolationSlice : public chi_mesh::FieldFunctionInterpolation
{
public:
  chi_mesh::Normal normal;
  chi_mesh::Normal binorm;
  chi_mesh::Normal tangent;
  chi_mesh::Vector3 point;

private:
  std::vector<int>                    intersecting_cell_indices;
  std::vector<FFICellIntersection*> cell_intersections;
  std::vector<int>                    cfem_local_nodes_needed_unmapped;
  std::vector<int>                    pwld_local_nodes_needed_unmapped;
  std::vector<int>                    pwld_local_cells_needed_unmapped;
public:
  FieldFunctionInterpolationSlice()
  {
    normal = chi_mesh::Normal(0.0,0.0,1.0);
    binorm = chi_mesh::Normal(0.0,1.0,0.0);
    tangent= chi_mesh::Normal(1.0,0.0,0.0);
  }

  //01
  void Initialize();

  //02
  void Execute();
private:
  void CFEMInterpolate(Vec field, std::vector<int> &mapping);
  void PWLDInterpolate(std::vector<double>& field, std::vector<int> &mapping);
public:
  //03
  void ExportPython(std::string base_name);
};


#endif