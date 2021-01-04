#ifndef _chi_angleaggregation_h
#define _chi_angleaggregation_h

#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiMesh/SweepUtilities/AngleSetGroup/anglesetgroup.h"

#include "ChiMath/Quadratures/angular_quadrature_base.h"

#include <memory>

//###################################################################
/**Angle aggregation has to cater for running the 8 corners of a 3D
 * partitioning, the 4 corners of a 2D partitioning (the latter 2 both being
 * polar angle aggregation) as well as single angle aggregation.
 *
 * At the most fundamental level this manifests as a number of angle indices
 * that share a SPDS, however SPDS do not have to be unique which allows
 * the notion of polar angle sets.
 * For single angle aggregation each SPDS is associated
 * with a single angle index. For polar angle aggregation a single SPDS can
 * have multiple indices associated with the same azimuthal angle but
 * different polar angles. We call this manifestation a "AngleSet".
 *
 * The octant based separation is achieved via the notion of "AngleSetGroup"
 * which will group angle sets for each quadrant or octant
 * (depending on 2D or 3D).*/
class chi_mesh::sweep_management::AngleAggregation
{
public:
  std::vector<AngleSetGroup*>  angle_set_groups;
  std::vector<SweepBndry*>     sim_boundaries;
  int                          number_of_groups=0;
  int                          number_of_group_subsets=0;
  std::shared_ptr<chi_math::AngularQuadrature> quadrature=nullptr;

private:
  bool is_setup=false;
  std::pair<int,int> number_angular_unknowns;
  bool num_ang_unknowns_avail = false;

public:
  chi_mesh::MeshContinuum* grid = nullptr;

  void Setup(const std::vector<SweepBndry*>& in_sim_boundaries,
             int in_number_of_groups,
             int in_number_of_group_subsets,
             std::shared_ptr<chi_math::AngularQuadrature>& in_quadrature,
             chi_mesh::MeshContinuum* in_grid);

public:
  double GetDelayedPsiNorm();
  void   ZeroOutgoingDelayedPsi();
  void   ZeroIncomingDelayedPsi();

  void InitializeReflectingBCs();
  void ResetReflectingBCs();

  std::pair<int,int> GetNumberOfAngularUnknowns();
  void AssembleAngularUnknowns(int& index, double* x_ref);
  void DisassembleAngularUnknowns(int& index, const double* x_ref);

};


#endif
