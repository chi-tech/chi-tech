#ifndef LINEAR_BOLTZMANN_SOLVER_H
#define LINEAR_BOLTZMANN_SOLVER_H

#include "A_LBSSolver/lbs_solver.h"

namespace lbs
{
//################################################################### Class def
/**A neutral particle transport solver.*/
class SteadyStateSolver : public LBSSolver
{
 public:
  //00
  explicit SteadyStateSolver(const std::string& in_text_name);
  ~SteadyStateSolver() override;

  SteadyStateSolver (const SteadyStateSolver&) = delete;
  SteadyStateSolver& operator= (const SteadyStateSolver&) = delete;

  //01
  void Initialize() override;
protected:
//  //01a
//  virtual void PerformInputChecks();
//  //01b
//  void PrintSimHeader();
//public:
//  //01c
//  void InitMaterials();
//protected:
//  //01d
//  virtual void InitializeSpatialDiscretization();
//  void ComputeUnitIntegrals();
//  //01e
//  void InitializeGroupsets();
//  //01f
//  void ComputeNumberOfMoments();
//  //01g
//  virtual void InitializeParrays();
//  //01h
//  void InitializeBoundaries();
//public:
//  //01i
//  void InitializePointSources();
protected:
  //01j
  virtual void InitializeSolverSchemes();
  virtual void InitializeWGSSolvers();




public:
  //02
  void Execute() override;

protected:

  //03a
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

  //03d
  void InitWGDSA(LBSGroupset& groupset);
public:
  void ExecuteWGDSA(LBSGroupset& groupset,
                    const std::vector<double>& ref_phi_old,
                    std::vector<double>& ref_phi_new);
  void AssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                   const std::vector<double>& ref_phi_old,
                                   const std::vector<double>& ref_phi_new,
                                   std::vector<double>& delta_phi_local);

  void AssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                   const std::vector<double>& phi_in,
                                   std::vector<double>& delta_phi_local);

  void DisAssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                      const std::vector<double>& delta_phi_local,
                                      std::vector<double>& ref_phi_new);
protected:
  static void CleanUpWGDSA(LBSGroupset& groupset);

  //03e
  void InitTGDSA(LBSGroupset& groupset);
public:
  void ExecuteTGDSA(LBSGroupset& groupset,
                    const std::vector<double>& ref_phi_old,
                    std::vector<double>& ref_phi_new);
  void AssembleTGDSADeltaPhiVector(const LBSGroupset& groupset,
                                   const std::vector<double>& ref_phi_old,
                                   const std::vector<double>& ref_phi_new,
                                   std::vector<double>& delta_phi_local);
  void AssembleTGDSADeltaPhiVector(const LBSGroupset& groupset,
                                   const std::vector<double>& phi_in,
                                   std::vector<double>& delta_phi_local);
  void DisAssembleTGDSADeltaPhiVector(const LBSGroupset& groupset,
                                      const std::vector<double>& delta_phi_local,
                                      std::vector<double>& ref_phi_new);
protected:
  static void CleanUpTGDSA(LBSGroupset& groupset);
  //03f
  void ResetSweepOrderings(LBSGroupset& groupset);

//  //04 File IO
//  //04a
//public:
//  void WriteRestartData(std::string folder_name, std::string file_base);
//  void ReadRestartData(std::string folder_name, std::string file_base);

//public:
//  //04b
//  void WriteGroupsetAngularFluxes(const LBSGroupset& groupset,
//                                  const std::string& file_base);
//  void ReadGroupsetAngularFluxes(LBSGroupset& groupset,
//                                 const std::string& file_base);
//
//  //04c
//  std::vector<double> MakeSourceMomentsFromPhi();
//  void WriteFluxMoments(const std::string& file_base,
//                        const std::vector<double>& flux_moments);
//  void ReadFluxMoments(const std::string& file_base,
//                       std::vector<double>& flux_moments,
//                       bool single_file=false);

//  //05a
//  void UpdateFieldFunctions();

  //Iterative Operations
  void SetSource(LBSGroupset& groupset,
                 std::vector<double>& destination_q,
                 const std::vector<double>& phi,
                 SourceFlags source_flags);
protected:
  double ComputePiecewiseChange(LBSGroupset& groupset);
  virtual std::shared_ptr<SweepChunk> SetSweepChunk(LBSGroupset& groupset);
  double ComputeFissionProduction(const std::vector<double>& phi);
public:
  virtual double ComputeFissionRate(bool previous);
protected:
  //Iterative Methods
//  bool Krylov(LBSGroupset& groupset,
//              chi_mesh::sweep_management::SweepScheduler& sweep_scheduler,
//              SourceFlags lhs_src_scope,
//              SourceFlags rhs_src_scope,
//              const SetSourceFunction& set_source_function,
//              bool log_info = true);

  //Vector assembly
public:
  void SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset, Vec x,
                                         const std::vector<double>& y,
                                         bool with_delayed_psi) override;
  void SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset, Vec x,
                                         PhiSTLOption which_phi) override;
  void SetGSSTLvectorFromPrimarySTLvector(LBSGroupset& groupset,
                                          std::vector<double>& x,
                                          const std::vector<double>& y,
                                          bool with_delayed_psi) override;
  void SetPrimarySTLvectorFromGSPETScVec(LBSGroupset& groupset, Vec x_src,
                                         std::vector<double>& y,
                                         bool with_delayed_psi) override;
  void SetPrimarySTLvectorFromGSSTLvector(LBSGroupset& groupset,
                                          const std::vector<double>& x_src,
                                          std::vector<double>& y,
                                          bool with_delayed_psi) override;
  void GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                     const std::vector<double>& x_src,
                                     std::vector<double>& y,
                                     bool with_delayed_psi) override;
//  void SetGroupScopedPETScVecFromPrimarySTLvector(int first_group_id,
//                                                  int last_group_id, Vec x,
//                                                  const std::vector<double>& y);
//  void SetPrimarySTLvectorFromGroupScopedPETScVec(
//    int first_group_id,
//    int last_group_id, Vec x_src,
//    std::vector<double>& y);
public:
  //compute_balance
  void ZeroOutflowBalanceVars(LBSGroupset& groupset);
  void ComputeBalance();

protected:
  //precursors
  void ComputePrecursors();
};

}

#endif
