#ifndef CHI_FFINTERPOLATION_H
#define CHI_FFINTERPOLATION_H

#include "../MeshContinuum/chi_meshcontinuum.h"
#include "../chi_mesh.h"
#include "ChiPhysics/FieldFunction/fieldfunction.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"

#define FFI_SLICE 1
#define FFI_LINE 2
#define FFI_VOLUME 3

#define FFI_SDM_CFEM    1
#define FFI_SDM_PWLD_GM 2   //PWLD with groups and moments

#define OP_SUM 10
#define OP_AVG 11
#define OP_MAX 12
#define OP_SUM_LUA 13
#define OP_AVG_LUA 14
#define OP_MAX_LUA 15

//###################################################################
/** Base class for field-function interpolation objects.*/
class chi_mesh::FieldFunctionInterpolation
{
public:
  chi_mesh::MeshContinuumPtr grid_view;

  std::vector<chi_physics::FieldFunction*> field_functions;

public:
  FieldFunctionInterpolation() : grid_view(nullptr) {}

  virtual void Initialize(){};
  virtual void Execute(){};
};



#endif