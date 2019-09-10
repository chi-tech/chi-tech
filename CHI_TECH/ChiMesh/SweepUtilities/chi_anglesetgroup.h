#ifndef _chi_angleset_group_h
#define _chi_angleset_group_h

#include "chi_angleset.h"



//###################################################################
/**Manages the workstages of a single angleset group.*/
class chi_mesh::SweepManagement::AngleSetGroup
{
public:
  std::vector<chi_mesh::SweepManagement::AngleSet*> angle_sets;

private:
  int current_angle_set = 0;


public:

  //###################################################################
  bool AngleSetGroupAdvance(chi_mesh::SweepManagement::SweepChunk *sweep_chunk,
                            int anglesetgroup_number)
  {
    //====================================== Return finished if angle sets
    //                                       depleted
    if (current_angle_set >= angle_sets.size())
    {
      return FLAG_FINISHED;
    }

    //====================================== Check anglesets can clear buffers
    for (int as=0; as<(current_angle_set-1); as++)
      angle_sets[as]->EnsureClearedBuffers();

    //====================================== Execute current angleset
    bool completion_status = angle_sets[current_angle_set]->
      AngleSetAdvance(sweep_chunk, current_angle_set +
      anglesetgroup_number * angle_sets.size());

    //====================================== Check if angle set finished
    if (completion_status ==  FLAG_FINISHED)
    {
      current_angle_set++;
      if (current_angle_set >= angle_sets.size())
      {
        return FLAG_FINISHED;
      }
      return FLAG_NOT_FINISHED;
    }

    return FLAG_NOT_FINISHED;
  }



public:
  //#################################################################
  void ResetSweep()
  {
    current_angle_set = 0;

    for (int as=0; as<angle_sets.size(); as++)
      angle_sets[as]->ResetSweepBuffers();
  }

};

#endif