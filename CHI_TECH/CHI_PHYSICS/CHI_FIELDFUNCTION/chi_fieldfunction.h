#ifndef _chi_field_function_h
#define _chi_field_function_h

#include "../chi_physics_namespace.h"
#include "../../CHI_MESH/CHI_MESHCONTINUUM/chi_meshcontinuum.h"
#include "CHI_MATH/SpatialDiscretization/spatial_discretization.h"

#include <petscksp.h>

#define FF_SDM_CFEM 1
#define FF_SDM_PWLD 2   //PWLD with groups and moments
#define FF_SDM_FV   3   //Finite Volume

class chi_physics::FieldFunction
{
public:
  std::string              text_name;
  int                      id;
  int                      type;
  int                      num_grps, num_moms;
  int                      grp, mom;

  chi_mesh::MeshContinuum* grid;
  SpatialDiscretization*      spatial_discretization;
  Vec                      field_vector;
  std::vector<double>*     field_vector_local;

  std::vector<int>*        local_cell_dof_array_address;

  FieldFunction()
  {
    id = 0;
    type = FF_SDM_CFEM;
    num_grps = 1;
    num_moms = 1;
    grp = 0;
    mom = 0;

    local_cell_dof_array_address = nullptr;
  }

  //01
  void ExportToVTK(std::string base_name, std::string field_name);
  void ExportToVTKG(std::string base_name, std::string field_name);
  //01a
  void ExportToVTKFV(std::string base_name, std::string field_name);
  void ExportToVTKFVG(std::string base_name, std::string field_name);
  //01b
  void ExportToVTKPWLC(std::string base_name, std::string field_name);
  void ExportToVTKPWLCG(std::string base_name, std::string field_name);
  //01c
  void ExportToVTKPWLD(std::string base_name, std::string field_name);
  void ExportToVTKPWLDG(std::string base_name, std::string field_name);

  void WritePVTU(std::string base_filename, std::string field_name, int num_grps=0);
};


#endif