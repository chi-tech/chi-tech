#ifndef _chi_sweepscheduler_H
#define _chi_sweepscheduler_H

#include "chi_angleaggregation.h"
#include "chi_sweepchunk_base.h"

//Scheduler types
#define FIFO           0
#define DEPTH_OF_GRAPH 1

typedef chi_mesh::SweepManagement::AngleSetGroup TAngleSetGroup;
typedef chi_mesh::SweepManagement::AngleSet      TAngleSet;
typedef chi_mesh::SweepManagement::STDG          TGSPO;
typedef std::vector<TGSPO*>                      TLEVELED_GRAPH;

//###################################################################
class chi_mesh::SweepManagement::SweepScheduler
{
private:
  int                       scheduler_type;
  AngleAggregation*         angle_agg;
  SweepChunk*               sweep_chunk;

  struct RULE_VALUES
  {
    TAngleSet* angle_set;
    int        depth_of_graph;
    int        sign_of_omegax;
    int        sign_of_omegay;
    int        sign_of_omegaz;
    int        set_index;
  };
  std::vector<RULE_VALUES> rule_values;

public:
  SweepScheduler(int in_scheduler_type,
                 AngleAggregation* in_angle_agg);

  void Sweep(SweepChunk* in_sweep_chunk=NULL);

private:
  void ScheduleAlgoFIFO();

  //02
  void InitializeAlgoDOG();
  void ScheduleAlgoDOG();
};

#endif