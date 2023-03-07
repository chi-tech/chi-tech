#ifndef CHI_SWEEP_BOUNDARY_BASE_H
#define CHI_SWEEP_BOUNDARY_BASE_H

#include "ChiMesh/chi_mesh.h"

#include <vector>
#include <limits>

namespace chi_math
{
  class AngularQuadrature;
}

namespace chi_mesh::sweep_management
{

enum class BoundaryType
{
  INCIDENT_ISOTROPIC_HOMOGENOUS     = 1, ///< One value for all angles, homogenous in space
  REFLECTING                        = 2, ///< Reflecting boundary condition about a normal
  INCIDENT_ANISOTROPIC_HETEROGENOUS = 3  ///< Complex different for each angle and face node
};

//###################################################################
/**Base class for sweep related boundaries.*/
class BoundaryBase
{
private:
  const chi_mesh::sweep_management::BoundaryType type;
protected:
  std::vector<double>  zero_boundary_flux;
  size_t num_groups;

public:
  explicit BoundaryBase(BoundaryType bndry_type,
                        size_t in_num_groups) :
    type(bndry_type),
    num_groups(in_num_groups)
  {
    zero_boundary_flux.resize(num_groups,0.0);
  }

  virtual ~BoundaryBase() = default;
  BoundaryType Type() const {return type;}
  bool     IsReflecting() const
  { return type == BoundaryType::REFLECTING; }


  virtual double* HeterogenousPsiIncoming(uint64_t cell_local_id,
                                          int face_num,
                                          int fi,
                                          int angle_num,
                                          int group_num,
                                          int gs_ss_begin);
  virtual double* HeterogenousPsiOutgoing(uint64_t cell_local_id,
                                          int face_num,
                                          int fi,
                                          int angle_num,
                                          int gs_ss_begin);
  virtual void UpdateAnglesReadyStatus(const std::vector<size_t>& angles,
                                       size_t gs_ss)
  {}
  virtual bool CheckAnglesReadyStatus(const std::vector<size_t>& angles,
                                      size_t gs_ss)
  {return true;}
  virtual void Setup(const chi_mesh::MeshContinuum& grid,
                     const chi_math::AngularQuadrature& quadrature) {}

  double* ZeroFlux(int group_num) {return &zero_boundary_flux[group_num];}
};

//###################################################################
/** Specified incident fluxes homogenous on a boundary.*/
class BoundaryIsotropicHomogenous : public BoundaryBase
{
private:
  std::vector<double> boundary_flux;
public:
  explicit
  BoundaryIsotropicHomogenous(size_t in_num_groups,
                              std::vector<double> ref_boundary_flux) :
    BoundaryBase(BoundaryType::INCIDENT_ISOTROPIC_HOMOGENOUS, in_num_groups),
    boundary_flux(std::move(ref_boundary_flux))
  {}

  double* HeterogenousPsiIncoming(
    uint64_t cell_local_id,
    int face_num,
    int fi,
    int angle_num,
    int group_num,
    int gs_ss_begin) override;
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
  BoundaryReflecting(size_t in_num_groups,
                     const chi_mesh::Normal& in_normal) :
    BoundaryBase(BoundaryType::REFLECTING,in_num_groups),
    normal(in_normal)
  {}

  double* HeterogenousPsiIncoming(uint64_t cell_local_id,
                                  int face_num,
                                  int fi,
                                  int angle_num,
                                  int group_num,
                                  int gs_ss_begin) override;
  double* HeterogenousPsiOutgoing(uint64_t cell_local_id,
                                  int face_num,
                                  int fi,
                                  int angle_num,
                                  int gs_ss_begin) override;

  void UpdateAnglesReadyStatus(const std::vector<size_t>& angles,
                               size_t gs_ss) override;
  bool CheckAnglesReadyStatus(const std::vector<size_t>& angles,
                              size_t gs_ss) override;
  void ResetAnglesReadyStatus();
};

class BoundaryFunction
{
public:
  virtual std::vector<double> Evaluate(
    size_t cell_global_id,
    int    cell_material_id,
    unsigned int face_index,
    unsigned int face_node_index,
    const chi_mesh::Vector3& face_node_location,
    const chi_mesh::Vector3& face_node_normal,
    const std::vector<int>& quadrature_angle_indices,
    const std::vector<chi_mesh::Vector3>& quadrature_angle_vectors,
    const std::vector<std::pair<double,double>>& quadrature_phi_theta_angles,
    const std::vector<int>& group_indices)
  {
    size_t num_angles = quadrature_angle_indices.size();
    size_t num_groups = group_indices.size();

    std::vector<double> psi(num_angles * num_groups, 0.0);

    return psi;
  }
};

//###################################################################
/** Specified incident fluxes homogenous on a boundary.*/
class BoundaryIncidentHeterogenous : public BoundaryBase
{
private:
  std::unique_ptr<BoundaryFunction> boundary_function;
  const uint64_t ref_boundary_id;

  typedef std::vector<double>       FaceNodeData;
  typedef std::vector<FaceNodeData> FaceData;
  typedef std::vector<FaceData>     CellData;

  std::vector<CellData> local_cell_data;
public:
  explicit
  BoundaryIncidentHeterogenous(size_t in_num_groups,
                               std::unique_ptr<BoundaryFunction> in_bndry_function,
                               uint64_t in_ref_boundary_id) :
    BoundaryBase(BoundaryType::INCIDENT_ANISOTROPIC_HETEROGENOUS, in_num_groups),
    boundary_function(std::move(in_bndry_function)),
    ref_boundary_id(in_ref_boundary_id)
  {}

  double* HeterogenousPsiIncoming(uint64_t cell_local_id,
                                  int face_num,
                                  int fi,
                                  int angle_num,
                                  int group_num,
                                  int gs_ss_begin) override;

  void Setup(const chi_mesh::MeshContinuum &grid,
             const chi_math::AngularQuadrature &quadrature) override;
};

}//namespace sweep_management

#endif //CHI_SWEEP_BOUNDARY_BASE_H