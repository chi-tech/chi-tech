#ifndef LINEAR_BOLTZMANN_SOLVER_H
#define LINEAR_BOLTZMANN_SOLVER_H

#include "A_LBSSolver/lbs_solver.h"

namespace lbs
{
//################################################################### Class def
/**A neutral particle transport solver.*/
class DiscOrdSteadyStateSolver : public LBSSolver
{
protected:
  typedef chi_mesh::sweep_management::SweepChunk SweepChunk;
protected:
  typedef std::shared_ptr<chi_math::AngularQuadrature> AngQuadPtr;
  typedef std::shared_ptr<chi_mesh::sweep_management::SPDS> SPDS_ptr;
  typedef std::vector<SPDS_ptr> SPDS_ptrs;
  typedef std::pair<UniqueSOGroupings, DirIDToSOMap> SwpOrderGroupingInfo;
  typedef chi_mesh::sweep_management::PRIMARY_FLUDS FLUDSTemplate;
  typedef std::shared_ptr<FLUDSTemplate> FLUDSTemplatePtr;
  typedef std::vector<FLUDSTemplatePtr> FLUDSTemplatePtrs;

  std::map<AngQuadPtr, SwpOrderGroupingInfo> quadrature_unq_so_grouping_map_;
  std::map<AngQuadPtr, SPDS_ptrs> quadrature_spds_map_;
  std::map<AngQuadPtr, FLUDSTemplatePtrs> quadrature_fluds_templates_map_;


 public:
  //00
  explicit DiscOrdSteadyStateSolver(const std::string& in_text_name);
  ~DiscOrdSteadyStateSolver() override;

  DiscOrdSteadyStateSolver (const DiscOrdSteadyStateSolver&) = delete;
  DiscOrdSteadyStateSolver& operator= (const DiscOrdSteadyStateSolver&) = delete;

  //01
  void Initialize() override;

protected:
  //01j
  void InitializeWGSSolvers() override;

public:
  //02
  void Execute() override;

  //03a
protected:
  void InitializeSweepDataStructures();
  //03aa
  static
  std::pair<UniqueSOGroupings, DirIDToSOMap>
  AssociateSOsAndDirections(const chi_mesh::MeshContinuum& grid,
                                 const chi_math::AngularQuadrature& quadrature,
                                 AngleAggregationType agg_type,
                                 lbs::GeometryType lbs_geo_type);
  //03b
  void InitFluxDataStructures(LBSGroupset& groupset);


protected:
  //03f
  void ResetSweepOrderings(LBSGroupset& groupset);

public:
  virtual std::shared_ptr<SweepChunk> SetSweepChunk(LBSGroupset& groupset);

protected:
  //Iterative Methods
  //Vector assembly
public:
  void SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset, Vec x,
                                         PhiSTLOption which_phi) override;

  void SetPrimarySTLvectorFromGSPETScVec(LBSGroupset& groupset, Vec x_src,
                                         PhiSTLOption which_phi) override;

  void GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                     PhiSTLOption from_which_phi,
                                     PhiSTLOption to_which_phi) override;

  void SetMultiGSPETScVecFromPrimarySTLvector(std::vector<int>& gs_ids,
                                              Vec x,
                                              PhiSTLOption which_phi) override;

  void SetPrimarySTLvectorFromMultiGSPETScVecFrom(std::vector<int>& gs_ids,
                                                  Vec x_src,
                                                  PhiSTLOption which_phi) override;

public:
  //compute_balance
  void ZeroOutflowBalanceVars(LBSGroupset& groupset);
  void ComputeBalance();
  std::vector<double>
  ComputeLeakage(int groupset_id, uint64_t boundary_id) const;

};

}

#endif
