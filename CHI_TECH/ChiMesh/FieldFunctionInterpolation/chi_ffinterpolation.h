#ifndef _chi_ffinterpolation_h
#define _chi_ffinterpolation_h

#include "../MeshContinuum/chi_meshcontinuum.h"
#include "../chi_mesh.h"
#include "ChiPhysics/FieldFunction/fieldfunction.h"

#define FFI_SLICE 1
#define FFI_LINE 2
#define FFI_VOLUME 3

#define FFI_SDM_CFEM    1
#define FFI_SDM_PWLD_GM 2   //PWLD with groups and moments

#define OP_SUM 10
#define OP_AVG 11
#define OP_MAX 12

//###################################################################
/** Base class for field-function interpolation objects.*/
class chi_mesh::FieldFunctionInterpolation
{
public:
  chi_mesh::MeshContinuum* grid_view;

  std::vector<chi_physics::FieldFunction*> field_functions;

  int spatial_discr_method;




public:
  FieldFunctionInterpolation()
  {
    grid_view = nullptr;
    spatial_discr_method = FFI_SDM_CFEM;
  }

  //01 utils
  bool CheckPlaneTetIntersect(chi_mesh::Normal plane_normal,
                              chi_mesh::Vector plane_point,
                              std::vector<chi_mesh::Vector*>* tet_points);
  bool CheckPlaneLineIntersect(chi_mesh::Normal plane_normal,
                               chi_mesh::Vector plane_point,
                               chi_mesh::Vector line_point_0,
                               chi_mesh::Vector line_point_1,
                               chi_mesh::Vector& intersection_point,
                               std::pair<double,double>& weights);
  bool CheckLineTriangleIntersect(std::vector<chi_mesh::Vector>& triangle_points,
                                  chi_mesh::Vector line_point_i,
                                  chi_mesh::Vector line_point_f);
  void CreateCFEMMapping(int num_grps, int num_moms, int g, int m,
                         Vec& x, Vec& x_cell,
                         std::vector<int>& cfem_nodes,
                         std::vector<int>* mapping);
  void CreatePWLDMapping(int num_grps, int num_moms, int g, int m,
                         std::vector<int>& pwld_nodes,
                         std::vector<int>& pwld_cells,
                         std::vector<int>& local_cell_dof_array_address,
                         std::vector<int>* mapping);
  void CreateFVMapping(int num_grps, int num_moms, int g, int m,
                       std::vector<int>& cells,
                       std::vector<int>* mapping);

//  void CreateCFEMMapping(int num_grps, int num_moms, int g, int m,
//                         Vec& x, Vec& x_cell,
//                         std::vector<int>& cfem_nodes,
//                         std::vector<int>* mapping);
  void CreatePWLDMapping(chi_physics::FieldFunction* field_function,
                         std::vector<int>& pwld_nodes,
                         std::vector<int>& pwld_cells,
                         std::vector<int>* mapping,int m=0,int g=0);
//  void CreateFVMapping(int num_grps, int num_moms, int g, int m,
//                       std::vector<int>& cells,
//                       std::vector<int>* mapping);

  virtual void Initialize(){};
  virtual void Execute(){};
};



#endif