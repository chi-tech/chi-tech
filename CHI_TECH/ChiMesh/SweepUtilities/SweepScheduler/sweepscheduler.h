#ifndef _chi_sweepscheduler_H
#define _chi_sweepscheduler_H

#include "ChiMesh/SweepUtilities/AngleAggregation/angleaggregation.h"
#include "ChiMesh/SweepUtilities/sweepchunk_base.h"

namespace chi_mesh::sweep_management
{
  enum class SchedulingAlgorithm {
    FIRST_IN_FIRST_OUT = 1,
    DEPTH_OF_GRAPH = 2
  };
}

typedef chi_mesh::sweep_management::AngleSetGroup TAngleSetGroup;
typedef chi_mesh::sweep_management::AngleSet      TAngleSet;
typedef chi_mesh::sweep_management::STDG          TGSPO;
typedef std::vector<TGSPO*>                      TLEVELED_GRAPH;

//###################################################################
class chi_mesh::sweep_management::SweepScheduler
{
private:
  SchedulingAlgorithm       scheduler_type;
  AngleAggregation*         angle_agg;
  SweepChunk*               sweep_chunk;
  const size_t              sweep_event_tag;

  struct RULE_VALUES
  {
    TAngleSet* angle_set;
    int        depth_of_graph;
    int        sign_of_omegax;
    int        sign_of_omegay;
    int        sign_of_omegaz;
    size_t     set_index;

    explicit RULE_VALUES(TAngleSet* ref_as) :
      angle_set(ref_as)
    {
      depth_of_graph = 0;
      set_index      = 0;
      sign_of_omegax = 1;
      sign_of_omegay = 1;
      sign_of_omegaz = 1;
    }
  };
  std::vector<RULE_VALUES> rule_values;

public:
  SweepScheduler(SchedulingAlgorithm in_scheduler_type,
                 AngleAggregation* in_angle_agg);

  void Sweep(SweepChunk* in_sweep_chunk=NULL);
  double GetAverageSweepTime();

private:
  void ScheduleAlgoFIFO();

  //02
  void InitializeAlgoDOG();
  void ScheduleAlgoDOG();
};

#endif