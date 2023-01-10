#ifndef CHI_SWEEP_BOUNDARY_BASE_H
#define CHI_SWEEP_BOUNDARY_BASE_H

#include "ChiMesh/chi_mesh.h"

#include <vector>
#include <limits>

namespace chi_mesh::sweep_management
{

enum class BoundaryType
{
  VACUUM                = 0,
  INCIDENT_HOMOGENOUS   = 1,
  REFLECTING            = 2,
  INCIDENT_HETEROGENOUS = 3
};

//###################################################################
/**Base class for sweep related boundaries.*/
class BoundaryBase
{
public:
  const chi_mesh::sweep_management::BoundaryType type;
  std::vector<double>& boundary_flux;
  std::vector<double>  zero_boundary_flux;

public:
  explicit BoundaryBase(BoundaryType bndry_type,
                        std::vector<double>& ref_boundary_flux) :
                        type(bndry_type),
                        boundary_flux(ref_boundary_flux)
  {
    zero_boundary_flux.resize(ref_boundary_flux.size(),0.0);
  }

  virtual ~BoundaryBase() = default;
  BoundaryType Type() const {return type;}
  bool     IsReflecting() const
  { return this->Type() == BoundaryType::REFLECTING; }


  virtual double* HeterogenousPsiIncoming(
                                  int angle_num,
                                  uint64_t cell_local_id,
                                  int face_num,
                                  int fi,
                                  int gs_ss_begin);
  virtual double* HeterogenousPsiOutgoing(
                                  int angle_num,
                                  uint64_t cell_local_id,
                                  int face_num,
                                  int fi,
                                  int gs_ss_begin);
  virtual void UpdateAnglesReadyStatus(const std::vector<size_t>& angles,
                                       size_t gs_ss)
  {}
  virtual bool CheckAnglesReadyStatus(const std::vector<size_t>& angles,
                                      size_t gs_ss)
  {return true;}

};

//###################################################################
/** Vacuum boundary.*/
class BoundaryVacuum : public BoundaryBase
{
public:
  explicit
  BoundaryVacuum(std::vector<double>& ref_boundary_flux) :
  BoundaryBase(BoundaryType::VACUUM,ref_boundary_flux)
  {}
};

//###################################################################
/** Specified incident fluxes homogenous on a boundary.*/
class BoundaryIncidentHomogenous : public BoundaryBase
{
public:
  explicit
  BoundaryIncidentHomogenous(std::vector<double>& ref_boundary_flux) :
    BoundaryBase(BoundaryType::INCIDENT_HOMOGENOUS,ref_boundary_flux)
  {}
};

//###################################################################
/** Specified incident fluxes homogenous on a boundary.*/
class BoundaryIncidentHeterogenous : public BoundaryBase
{
public:
  explicit
  BoundaryIncidentHeterogenous(std::vector<double>& ref_boundary_flux) :
    BoundaryBase(BoundaryType::INCIDENT_HETEROGENOUS,ref_boundary_flux)
  {}
};

//###################################################################
/** Reflective boundary condition.*/
class BoundaryReflecting : public BoundaryBase
{
public:
  const chi_mesh::Normal normal;
  bool  opposing_reflected = false;

  typedef std::vector<double> DOFVec;   //Groups per DOF
  typedef std::vector<DOFVec> FaceVec;  //DOFs per face
  typedef std::vector<FaceVec> CellVec; //Faces per cell
  typedef std::vector<CellVec> AngVec;  //Cell per angle

  //angle,cell,face,dof,group
  //Populated by angle aggregation
  std::vector<AngVec>              hetero_boundary_flux;
  std::vector<AngVec>              hetero_boundary_flux_old;
  double                           pw_change=0.0;

  std::vector<int>                 reflected_anglenum;
  std::vector<std::vector<bool>>   angle_readyflags;

public:
  BoundaryReflecting(std::vector<double>& ref_boundary_flux,
                     const chi_mesh::Normal& in_normal) :
  BoundaryBase(BoundaryType::REFLECTING,ref_boundary_flux),
  normal(in_normal)
  {}

  double* HeterogenousPsiIncoming(
                          int angle_num,
                          uint64_t cell_local_id,
                          int face_num,
                          int fi,
                          int gs_ss_begin) override;
  double* HeterogenousPsiOutgoing(
                          int angle_num,
                          uint64_t cell_local_id,
                          int face_num,
                          int fi,
                          int gs_ss_begin) override;

  void UpdateAnglesReadyStatus(const std::vector<size_t>& angles,
                               size_t gs_ss) override;
  bool CheckAnglesReadyStatus(const std::vector<size_t>& angles,
                              size_t gs_ss) override;
  void ResetAnglesReadyStatus();
};
}//namespace sweep_management

#endif //CHI_SWEEP_BOUNDARY_BASE_H