#include "sweepscheduler.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

#include <sstream>

//###################################################################
/**Initializes the Depth-Of-Graph algorithm.*/
void chi_mesh::sweep_management::SweepScheduler::InitializeAlgoDOG()
{
  //================================================== Load all anglesets
  //                                                   in preperation for
  //                                                   sorting
  //======================================== Loop over angleset groups
  for (size_t q=0; q<angle_agg->angle_set_groups.size(); q++)
  {
    TAngleSetGroup* angleset_group = angle_agg->angle_set_groups[q];

    //================================= Loop over anglesets in group
    size_t num_anglesets = angleset_group->angle_sets.size();
    for (size_t as=0; as<num_anglesets; as++)
    {
      TAngleSet* angleset           = angleset_group->angle_sets[as];
      auto       spds               = angleset->GetSPDS();
      TLEVELED_GRAPH& leveled_graph = spds->global_sweep_planes;

      //========================== Find location depth
      size_t loc_depth = -1;
      for (size_t level=0; level<leveled_graph.size(); level++)
      {
        for (size_t index=0; index<leveled_graph[level]->item_id.size(); index++)
        {
          if (leveled_graph[level]->item_id[index] == chi_mpi.location_id)
          {
            loc_depth = leveled_graph.size()-level;
            break;
          }
        }//for locations in plane
      }//for sweep planes

      //========================== Set up rule values
      if (loc_depth>=0)
      {
        RULE_VALUES new_rule_vals(angleset);
        new_rule_vals.depth_of_graph = loc_depth;
        new_rule_vals.set_index      = as + q * num_anglesets;

        new_rule_vals.sign_of_omegax = (spds->omega.x >= 0)?2:1;
        new_rule_vals.sign_of_omegay = (spds->omega.y >= 0)?2:1;
        new_rule_vals.sign_of_omegaz = (spds->omega.z >= 0)?2:1;

        rule_values.push_back(new_rule_vals);
      }
      else
      {
        chi_log.Log(LOG_ALLERROR)
          << "Location depth not found in Depth-Of-Graph algorithm.";
        exit(EXIT_FAILURE);
      }

    }//for anglesets
  }//for quadrants/anglesetgroups

  //================================================== Init sort functions
  struct
  {
    bool operator()(RULE_VALUES a, RULE_VALUES b)
    {
      return a.depth_of_graph > b.depth_of_graph;
    }
  }compare_D;

  struct
  {
    bool operator()(RULE_VALUES a, RULE_VALUES b)
    {
      return (a.depth_of_graph == b.depth_of_graph) and
             (a.sign_of_omegax > b.sign_of_omegax);
    }
  }compare_omega_x;

  struct
  {
    bool operator()(RULE_VALUES a, RULE_VALUES b)
    {
      return (a.depth_of_graph == b.depth_of_graph) and
             (a.sign_of_omegax == b.sign_of_omegax) and
             (a.sign_of_omegay > b.sign_of_omegay);
    }
  }compare_omega_y;

  struct
  {
    bool operator()(RULE_VALUES a, RULE_VALUES b)
    {
      return (a.depth_of_graph == b.depth_of_graph) and
             (a.sign_of_omegax == b.sign_of_omegax) and
             (a.sign_of_omegay == b.sign_of_omegay) and
             (a.sign_of_omegaz > b.sign_of_omegaz);
    }
  }compare_omega_z;

  //================================================== Sort
  std::stable_sort(rule_values.begin(),rule_values.end(),compare_D);
  std::stable_sort(rule_values.begin(),rule_values.end(),compare_omega_x);
  std::stable_sort(rule_values.begin(),rule_values.end(),compare_omega_y);
  std::stable_sort(rule_values.begin(),rule_values.end(),compare_omega_z);

}

//###################################################################
/**Executes the Depth-Of-Graph algorithm.*/
void chi_mesh::sweep_management::SweepScheduler::ScheduleAlgoDOG()
{
  typedef ExecutionPermission ExePerm;
  typedef AngleSetStatus Status;

  chi_log.LogEvent(sweep_event_tag, ChiLog::EventType::EVENT_BEGIN);

  auto ev_info =
    std::make_shared<ChiLog::EventInfo>(std::string("Sweep initiated"));

  chi_log.LogEvent(sweep_event_tag,
                   ChiLog::EventType::SINGLE_OCCURRENCE,ev_info);

  //==================================================== Loop till done
  bool finished = false;
  size_t scheduled_angleset = 0;
  while (!finished)
  {
    finished = true;
    for (size_t as=0; as<rule_values.size(); as++)
    {
      TAngleSet* angleset = rule_values[as].angle_set;
      int angset_number = rule_values[as].set_index;

      //=============================== Query angleset status
      // Status will here be one of the following:
      //  - RECEIVING.
      //      Meaning it is either waiting for messages or actively receiving it
      //  - READY_TO_EXECUTE.
      //      Meaning it has received all upstream data and can be executed
      //  - FINISHED.
      //      Meaning the angleset has executed its sweep chunk
      Status status = angleset->
        AngleSetAdvance(sweep_chunk,
                        angset_number,
                        sweep_timing_events_tag,
                        ExePerm::NO_EXEC_IF_READY);

      //=============================== Execute if ready and allowed
      // If this angleset is the one scheduled to run
      // and it is ready then it will be given permission
      if (status == Status::READY_TO_EXECUTE /*and as == scheduled_angleset*/)
      {
        std::stringstream message_i;
        message_i
          << "Angleset " << angset_number
          << " executed on location " << chi_mpi.location_id;

        auto ev_info_i = std::make_shared<ChiLog::EventInfo>(message_i.str());

        chi_log.LogEvent(sweep_event_tag,
                         ChiLog::EventType::SINGLE_OCCURRENCE,ev_info_i);

        status = angleset->
          AngleSetAdvance(sweep_chunk,
                          angset_number,
                          sweep_timing_events_tag,
                          ExePerm::EXECUTE);

        std::stringstream message_f;
        message_f
          << "Angleset " << angset_number
          << " finished on location " << chi_mpi.location_id;

        auto ev_info_f = std::make_shared<ChiLog::EventInfo>(message_f.str());

        chi_log.LogEvent(sweep_event_tag,
                         ChiLog::EventType::SINGLE_OCCURRENCE,ev_info_f);

        scheduled_angleset++; //Schedule the next angleset
      }

      if (status != Status::FINISHED)
        finished = false;
    }//for each angleset rule
  }//while not finished

  //================================================== Reset all
  for (auto angset_group : angle_agg->angle_set_groups)
    angset_group->ResetSweep();

  for (auto bndry : angle_agg->sim_boundaries)
  {
    if (bndry->Type() == chi_mesh::sweep_management::BoundaryType::REFLECTING)
    {
      auto rbndry = (chi_mesh::sweep_management::BoundaryReflecting*)bndry;
      rbndry->ResetAnglesReadyStatus();
    }
  }

  //================================================== Receive delayed data
  MPI_Barrier(MPI_COMM_WORLD);
  for (auto sorted_angleset : rule_values)
  {
    TAngleSet *angleset = sorted_angleset.angle_set;
    angleset->ReceiveDelayedData(sorted_angleset.set_index);
  }

  chi_log.LogEvent(sweep_event_tag, ChiLog::EventType::EVENT_END);
}