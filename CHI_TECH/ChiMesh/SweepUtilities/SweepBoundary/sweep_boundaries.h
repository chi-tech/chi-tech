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

public:
  explicit BoundaryBase(BoundaryType bndry_type,
                        std::vector<double>& ref_boundary_flux) :
                        type(bndry_type),
                        boundary_flux(ref_boundary_flux)
  { }

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

}

#endif