#ifndef SPATIAL_DISCRETIZATION_H
#define SPATIAL_DISCRETIZATION_H

#include "ChiMesh/chi_mesh.h"
#include "../Quadratures/quadrature.h"
#include "ChiMath/chi_math.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <petscksp.h>

#define DEG1         0
#define DEG3         1
#define DEG3_SURFACE 2

class SpatialDiscretization;

typedef std::shared_ptr<SpatialDiscretization> SpatialDiscretizationPtr;

class SpatialDiscretization
{
public:
  int dim;
  const chi_math::SpatialDiscretizationType type;

  chi_mesh::MeshContinuumPtr ref_grid;

protected:
  typedef chi_math::SpatialDiscretizationType SDMType;

protected:
  //00
  explicit SpatialDiscretization(int dim,
                                 chi_mesh::MeshContinuumPtr in_grid,
                                 SDMType in_type =
                                          SDMType::UNDEFINED);

protected:
  //01
  virtual void PreComputeCellSDValues(chi_mesh::MeshContinuumPtr grid);

public:
  virtual
  void BuildSparsityPattern(chi_mesh::MeshContinuumPtr grid,
                            std::vector<int>& nodal_nnz_in_diag,
                            std::vector<int>& nodal_nnz_off_diag,
                            chi_math::UnknownManager& unknown_manager)
                            {}

  virtual
  size_t GetNumLocalDOFs(chi_mesh::MeshContinuumPtr grid,
                         chi_math::UnknownManager& unknown_manager)
                         {return 0;}
  virtual
  size_t GetNumGlobalDOFs(chi_mesh::MeshContinuumPtr grid,
                          chi_math::UnknownManager& unknown_manager)
                          {return 0;}

protected:
  //02
  /**Develops a localized view of a petsc vector.
   * Each spatial discretization has a specialization of this
   * method.*/
  virtual void LocalizePETScVector(Vec petsc_vector,
                                   std::vector<double>& local_vector,
                                   chi_math::UnknownManager& unknown_manager)
  {}

};

#endif