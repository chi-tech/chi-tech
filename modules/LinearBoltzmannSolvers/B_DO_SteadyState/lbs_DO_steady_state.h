#ifndef LINEAR_BOLTZMANN_SOLVER_H
#define LINEAR_BOLTZMANN_SOLVER_H

#include "A_LBSSolver/lbs_solver.h"

namespace lbs
{
//################################################################### Class def
/**A neutral particle transport solver.*/
class DiscOrdSteadyStateSolver : public LBSSolver
{
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
  void ComputeSweepOrderings(LBSGroupset& groupset) const;
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

protected:
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

public:
  //compute_balance
  void ZeroOutflowBalanceVars(LBSGroupset& groupset);
  void ComputeBalance();

};

}

#endif
