#include "chi_sweepscheduler.h"

#include <chi_mpi.h>
#include <chi_log.h>
#include <ChiTimer/chi_timer.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

extern double chi_global_timings[20];



//###################################################################
/**Initializes the Depth-Of-Graph algorithm.*/
void chi_mesh::SweepManagement::SweepScheduler::InitializeAlgoDOG()
{
  //================================================== Load all anglesets
  //                                                   in preperation for
  //                                                   sorting
  //======================================== Loop over angleset groups
  for (int q=0; q<angle_agg->angle_set_groups.size(); q++)
  {
    TAngleSetGroup* angleset_group = angle_agg->angle_set_groups[q];

    //================================= Loop over anglesets in group
    for (int as=0; as<angleset_group->angle_sets.size(); as++)
    {
      TAngleSet* angleset           = angleset_group->angle_sets[as];
      TLEVELED_GRAPH& leveled_graph = angleset->GetSPDS()->global_sweep_planes;

      //========================== Find location depth
      int depth_counter=0;
      for (int level=0; level<leveled_graph.size(); level++)
      {
        for (int index=0; index<leveled_graph[level]->item_id.size(); index++)
        {
          depth_counter++;
        }//for locations in plane
      }//for sweep planes

      int max_depth = depth_counter;
      depth_counter=0;
      int loc_depth = -1;
      for (int level=0; level<leveled_graph.size(); level++)
      {
        for (int index=0; index<leveled_graph[level]->item_id.size(); index++)
        {
          depth_counter++;
          if (leveled_graph[level]->item_id[index] == chi_mpi.location_id)
          {
//            loc_depth = max_depth - depth_counter;
            loc_depth = leveled_graph.size()-level;
          }//if location
        }//for locations in plane
      }//for sweep planes

      if (loc_depth>=0)
      {
        RULE_VALUES new_rule_vals;
        new_rule_vals.angle_set      = angleset;
        new_rule_vals.depth_of_graph = loc_depth;
        new_rule_vals.set_index      = as +
                                       q * angleset_group->angle_sets.size();

        if (angleset->GetSPDS()->omega.x >= 0)
          new_rule_vals.sign_of_omegax = 2;
        else
          new_rule_vals.sign_of_omegax = 1;

        if (angleset->GetSPDS()->omega.y >= 0)
          new_rule_vals.sign_of_omegay = 2;
        else
          new_rule_vals.sign_of_omegay = 1;

        if (angleset->GetSPDS()->omega.z >= 0)
          new_rule_vals.sign_of_omegaz = 2;
        else
          new_rule_vals.sign_of_omegaz = 1;

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

  //==================================================== Init sort functions
  struct
  {
    bool operator()(RULE_VALUES a, RULE_VALUES b)
    {
      if (a.depth_of_graph > b.depth_of_graph)
      {return true;}
      else
      {return false;}
    }
  }compare_D;

  struct
  {
    bool operator()(RULE_VALUES a, RULE_VALUES b)
    {
      if ( (a.depth_of_graph == b.depth_of_graph) &&
           (a.sign_of_omegax > b.sign_of_omegax) )
      {return true;}
      else
      {return false;}
    }
  }compare_omega_x;

  struct
  {
    bool operator()(RULE_VALUES a, RULE_VALUES b)
    {
      if ( (a.depth_of_graph == b.depth_of_graph) &&
           (a.sign_of_omegax == b.sign_of_omegax) &&
           (a.sign_of_omegay  > b.sign_of_omegay) )
      {return true;}
      else
      {return false;}
    }
  }compare_omega_y;

  struct
  {
    bool operator()(RULE_VALUES a, RULE_VALUES b)
    {
      if ( (a.depth_of_graph == b.depth_of_graph) &&
           (a.sign_of_omegax == b.sign_of_omegax) &&
           (a.sign_of_omegay == b.sign_of_omegay) &&
           (a.sign_of_omegaz  > b.sign_of_omegaz)  )
      {return true;}
      else
      {return false;}
    }
  }compare_omega_z;

  //==================================================== Sort
  std::stable_sort(rule_values.begin(),rule_values.end(),compare_D);
  std::stable_sort(rule_values.begin(),rule_values.end(),compare_omega_x);
  std::stable_sort(rule_values.begin(),rule_values.end(),compare_omega_y);
  std::stable_sort(rule_values.begin(),rule_values.end(),compare_omega_z);

}

//###################################################################
/**Executes the Depth-Of-Graph algorithm.*/
void chi_mesh::SweepManagement::SweepScheduler::ScheduleAlgoDOG()
{
  ChiTimer t16_sweeptime; t16_sweeptime.Reset();

  //==================================================== Loop till done
  bool completion_status = FLAG_FINISHED;
  bool first_group = true;
  while ((completion_status == FLAG_NOT_FINISHED) || (first_group))
  {
    completion_status = FLAG_FINISHED;
    first_group=false;
    for (int as=0; as<rule_values.size(); as++)
    {
      TAngleSet* angleset = rule_values[as].angle_set;

      angleset->EnsureClearedBuffers();
      bool angleset_status = angleset->
        AngleSetAdvance(sweep_chunk, rule_values[as].set_index);

      if (angleset_status == FLAG_NOT_FINISHED)
      {
        completion_status = FLAG_NOT_FINISHED;
        //break;
      }

    }
  }

  //================================================== Reset all
  for (int q=0; q<angle_agg->angle_set_groups.size(); q++)
  {
    angle_agg->angle_set_groups[q]->ResetSweep();
  }

  MPI_Barrier(MPI_COMM_WORLD);
  for (int as=0; as<rule_values.size(); as++)
  {
    TAngleSet *angleset = rule_values[as].angle_set;
    angleset->ReceiveDelayedData(rule_values[as].set_index);
  }
//  chi_log.Log(LOG_ALL) << "Done with Sweep";
//  exit(0);
  chi_global_timings[16] += t16_sweeptime.GetTime()/1000.0;
  chi_global_timings[17] += 1.0;

}