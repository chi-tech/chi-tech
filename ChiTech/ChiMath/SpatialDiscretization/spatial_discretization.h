#ifndef SPATIAL_DISCRETIZATION_H
#define SPATIAL_DISCRETIZATION_H

#include "ChiMesh/chi_mesh.h"
#include "../Quadratures/quadrature.h"
#include "ChiMath/chi_math.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <petscksp.h>

class SpatialDiscretization
{
public:
  const int dim;
  const chi_math::SpatialDiscretizationType type;

  const chi_mesh::MeshContinuumPtr ref_grid;
  const chi_math::CoordinateSystemType cs_type;

protected:
  typedef chi_math::SpatialDiscretizationType SDMType;

protected:
  //00
  explicit
  SpatialDiscretization(int in_dim,
                        chi_mesh::MeshContinuumPtr& in_grid,
                        chi_math::CoordinateSystemType in_cs_type,
                        SDMType in_type = SDMType::UNDEFINED) :
    dim(in_dim),
    type(in_type),
    ref_grid(in_grid),
    cs_type(in_cs_type)
  {
  }

protected:
  //01
  virtual void PreComputeCellSDValues() {}

public:
  virtual
  void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                            std::vector<int64_t>& nodal_nnz_off_diag,
                            chi_math::UnknownManager& unknown_manager)
                            {}

  virtual
  int64_t MapDOF(const chi_mesh::Cell& cell,
                 unsigned int node,
                 const chi_math::UnknownManager& unknown_manager,
                 unsigned int unknown_id,
                 unsigned int component) const {return 0;}

  virtual
  int64_t MapDOFLocal(const chi_mesh::Cell& cell,
                      unsigned int node,
                      const chi_math::UnknownManager& unknown_manager,
                      unsigned int unknown_id,
                      unsigned int component) const {return 0;}

  virtual
  int64_t MapDOF(const chi_mesh::Cell& cell, unsigned int node) const {return 0;}
  virtual
  int64_t MapDOFLocal(const chi_mesh::Cell& cell, unsigned int node) const {return 0;}

  virtual
  size_t GetNumLocalDOFs(chi_math::UnknownManager& unknown_manager)
                         {return 0;}
  virtual
  size_t GetNumGlobalDOFs(chi_math::UnknownManager& unknown_manager)
                          {return 0;}

  virtual
  size_t GetCellNumNodes(const chi_mesh::Cell& cell) = 0;

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