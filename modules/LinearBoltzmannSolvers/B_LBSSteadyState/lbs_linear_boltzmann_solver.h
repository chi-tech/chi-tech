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


protected:
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
