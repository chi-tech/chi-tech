#ifndef _chi_angleaggregation_h
#define _chi_angleaggregation_h

#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiMesh/SweepUtilities/AngleSetGroup/anglesetgroup.h"

#include "ChiMath/Quadratures/product_quadrature.h"

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
  chi_math::ProductQuadrature* quadrature=nullptr;

public:
  chi_mesh::MeshContinuum* grid;

public:
  double GetDelayedPsiNorm();
  void   ResetDelayedPsi();

  void InitializeReflectingBCs();
  void ResetReflectingBCs();

};


#endif
