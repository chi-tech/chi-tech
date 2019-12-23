#ifndef _chi_sweep_bndry_base_h
#define _chi_sweep_bndry_base_h

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

  virtual ~BoundaryBase() {}
  const    BoundaryType Type() {return type;}
};

//###################################################################
/** Vacuum boundary.*/
class BoundaryVacuum : public BoundaryBase
{
public:
  BoundaryVacuum(std::vector<double>& ref_boundary_flux) :
  BoundaryBase(BoundaryType::VACUUM,ref_boundary_flux)
  {}
};

//###################################################################
/** Specified incident fluxes homogenous on a boundary.*/
class BoundaryIncidentHomogenous : public BoundaryBase
{
public:
  BoundaryIncidentHomogenous(std::vector<double>& ref_boundary_flux) :
    BoundaryBase(BoundaryType::INCIDENT_HOMOGENOUS,ref_boundary_flux)
  {}
};

//###################################################################
/** Specified incident fluxes homogenous on a boundary.*/
class BoundaryIncidentHeterogenous : public BoundaryBase
{
public:
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

  typedef std::vector<double> DOFVec;   //Groups per DOF
  typedef std::vector<DOFVec> FaceVec;  //DOFs per face
  typedef std::vector<FaceVec> CellVec; //Faces per cell
  typedef std::vector<CellVec> AngVec;  //Cell per angle

  //angle,cell,face,dof,group
  //Populated by angle aggregation
  std::vector<AngVec>              hetero_boundary_flux;
  std::vector<int>                 reflected_anglenum;
  std::vector<std::vector<bool>>   angle_readyflags;

public:
  BoundaryReflecting(std::vector<double>& ref_boundary_flux,
                     chi_mesh::Normal in_normal) :
  BoundaryBase(BoundaryType::REFLECTING,ref_boundary_flux),
  normal(in_normal)
  {}

  /**Sets flags indicating reflected angles are ready to execute.*/
  void UpdateAnglesReadyStatus(std::vector<int> angles, int gs_ss)
  {
    for (auto& n : angles)
      angle_readyflags[reflected_anglenum[n]][gs_ss] = true;
  }

  /**Checks to see if angles are ready to execute.*/
  bool CheckAnglesReadyStatus(std::vector<int> angles, int gs_ss)
  {
    bool ready_flag = true;
    for (auto& n : angles)
      if (hetero_boundary_flux[reflected_anglenum[n]].size()>0)
        if (not angle_readyflags[n][gs_ss]) return false;

    return ready_flag;
  }

  /**Resets angle ready flags to false.*/
  void ResetAnglesReadyStatus()
  {
    for (auto& flags : angle_readyflags)
      for (int gs_ss=0; gs_ss<flags.size(); ++gs_ss)
        flags[gs_ss] = false;
  }
};

}

#endif