#ifndef CHITECH_LBS_DISCRETE_ORDINATES_SOLVER_H
#define CHITECH_LBS_DISCRETE_ORDINATES_SOLVER_H

#include "A_LBSSolver/lbs_solver.h"

namespace lbs
{

/**Base class for Discrete Ordinates solvers. This class mostly establishes
 * utilities related to sweeping. From here we can derive a steady-state,
 * transient, adjoint, and k-eigenvalue solver.*/
class LBSDiscreteOrdinatesSolver : public LBSSolver
{
protected:
  typedef chi_mesh::sweep_management::SweepChunk SweepChunk;

  typedef std::shared_ptr<chi_math::AngularQuadrature> AngQuadPtr;

  typedef std::pair<UniqueSOGroupings, DirIDToSOMap> SwpOrderGroupingInfo;

  typedef std::shared_ptr<chi_mesh::sweep_management::SPDS> SPDS_ptr;
  typedef std::vector<SPDS_ptr> SPDS_ptrs;

  typedef chi_mesh::sweep_management::PRIMARY_FLUDS FLUDSTemplate;
  typedef std::shared_ptr<FLUDSTemplate> FLUDSTemplatePtr;
  typedef std::vector<FLUDSTemplatePtr> FLUDSTemplatePtrs;

protected:
  std::map<AngQuadPtr, SwpOrderGroupingInfo> quadrature_unq_so_grouping_map_;
  std::map<AngQuadPtr, SPDS_ptrs> quadrature_spds_map_;
  std::map<AngQuadPtr, FLUDSTemplatePtrs> quadrature_fluds_templates_map_;

public:
  static chi_objects::InputParameters GetInputParameters();
  explicit LBSDiscreteOrdinatesSolver(
    const chi_objects::InputParameters& params);
protected:
  explicit LBSDiscreteOrdinatesSolver(const std::string& text_name);

public:
  virtual ~LBSDiscreteOrdinatesSolver() override;

  std::pair<size_t, size_t> GetNumPhiIterativeUnknowns() override;

  // 01
  void Initialize() override;

protected:
  // 01j
  void InitializeWGSSolvers() override;

  // Sweep Data
  void InitializeSweepDataStructures();
  static std::pair<UniqueSOGroupings, DirIDToSOMap>
  AssociateSOsAndDirections(const chi_mesh::MeshContinuum& grid,
                            const chi_math::AngularQuadrature& quadrature,
                            AngleAggregationType agg_type,
                            lbs::GeometryType lbs_geo_type);
  void InitFluxDataStructures(LBSGroupset& groupset);
  void ResetSweepOrderings(LBSGroupset& groupset);
  virtual std::shared_ptr<SweepChunk> SetSweepChunk(LBSGroupset& groupset);

  // Vector assembly
public:
  void ScalePhiVector(PhiSTLOption which_phi, double value) override;
  void SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset,
                                         Vec x,
                                         PhiSTLOption which_phi) override;

  void SetPrimarySTLvectorFromGSPETScVec(LBSGroupset& groupset,
                                         Vec x_src,
                                         PhiSTLOption which_phi) override;

  void GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                     PhiSTLOption from_which_phi,
                                     PhiSTLOption to_which_phi) override;

  void SetMultiGSPETScVecFromPrimarySTLvector(const std::vector<int>& gs_ids,
                                              Vec x,
                                              PhiSTLOption which_phi) override;

  void SetPrimarySTLvectorFromMultiGSPETScVecFrom(
    const std::vector<int>& gs_ids, Vec x_src, PhiSTLOption which_phi) override;

  // compute_balance
public:
  void ZeroOutflowBalanceVars(LBSGroupset& groupset);
  void ComputeBalance();

  // compute leakage
public:
  std::vector<double> ComputeLeakage(int groupset_id,
                                     uint64_t boundary_id) const;
};

} // namespace lbs

#endif // CHITECH_LBS_DISCRETE_ORDINATES_SOLVER_H
