#ifndef _k_eigen_solver_h
#define _k_eigen_solver_h

#include "LinearBoltzmannSolver/lbs_linear_boltzmann_solver.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <string>

namespace  LinearBoltzmann { namespace KEigenvalue
{

/**A k-eigenvalue neutron transport solver.*/
class Solver : public LinearBoltzmann::Solver
{
private:
  size_t source_event_tag;

public:

  // Total number of precursors on all materials
  size_t num_precursors;

  // Current estimate of the k-eigenvalue
  double k_eff = 1.0;

  /** This is a local-to-global index mapping for the precursors
  //  defined on each material. The precursor vector Nj_new_local
  //  stores all precursor concentrations defined in the problem
  //  on each node, whether a precursor exists on the material a
  //  node is located on. This necessitates a global numbering
  //  system for the precursors. The TransportCrossSections
  // class defines all precursor information locally, with indices
  // from 0 to num_precursors - 1 for that material. The precursor_map
  // is a vector of vectors containing the global numbering of the
  // precursor species. The outer vector is of the same length as
  // the number of material_xs vector. Each of the inner vectors have
  // length corresponding to num_precursors for the corresponding
  // TransportCrossSection object.
  **/
  std::vector<std::vector<size_t>> precursor_map;

  // Additional phi vector
  std::vector<double> phi_prev_local;

  // Precursor vector and unknown manager
  std::vector<double> Nj_new_local;
  chi_math::UnknownManager Nj_unk_man;

  // IterativeMethods
  void PowerIteration();
  
  // Iterative operations
  void SetKSource(LBSGroupset& groupset,
                  std::vector<double>& destination_q,
                  SourceFlags source_flags);
  double ComputeProduction();
  void InitializePrecursors();

  // Execute method
  void InitializeKSolver();
  void ExecuteKSolver();
  
};

} }

#endif