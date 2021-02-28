#ifndef CHI_FFINTER_SLICE_H
#define CHI_FFINTER_SLICE_H

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
  int v0_g_index=0, v1_g_index=0;
  int v0_dofindex_cell=0, v1_dofindex_cell=0;
  chi_mesh::Vector3 point;
  chi_mesh::Vector3 point2d;
  std::pair<double,double> weights;
  double point_value=0.0;

  FFIFaceEdgeIntersection()
  {
    point_value = 0.0;
  }
};

struct FFICellIntersection
{
  int cell_local_index=0;
  std::vector<FFIFaceEdgeIntersection> intersections;
  chi_mesh::Vector3 intersection_centre;
  chi_mesh::Vector3 intersection_2d_centre;
  double cell_avg_value=0.0;

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
  chi_mesh::Normal normal =chi_mesh::Normal(0.0,0.0,1.0);
  chi_mesh::Normal binorm =chi_mesh::Normal(0.0,1.0,0.0);
  chi_mesh::Normal tangent=chi_mesh::Normal(1.0,0.0,0.0);
  chi_mesh::Vector3 point;

private:
  std::vector<uint64_t>               intersecting_cell_indices;
  std::vector<FFICellIntersection>    cell_intersections;
  std::vector<int>                    cfem_local_nodes_needed_unmapped;
  std::vector<int>                    cfem_local_cells_needed_unmapped;
  std::vector<int>                    pwld_local_nodes_needed_unmapped;
  std::vector<int>                    pwld_local_cells_needed_unmapped;
public:
  FieldFunctionInterpolationSlice() = default;

  //01
  void Initialize() override;

  //02
  void Execute() override;
private:
  void CFEMInterpolate(Vec field, std::vector<uint64_t> &mapping);
  void PWLDInterpolate(std::vector<double>& field, std::vector<uint64_t>& mapping);
public:
  //03
  void ExportPython(std::string base_name);
};


#endif