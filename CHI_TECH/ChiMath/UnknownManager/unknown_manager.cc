#include "unknown_manager.h"

#include "chi_log.h"

//###################################################################
/**Adds an unknown to the manager. This method will figure out
 * where the last unknown ends and where to begin the next one.*/
unsigned int chi_math::UnknownManager::
  AddUnknown(const UnknownType unk_type,
             const unsigned int dimension)
{
  auto& log = ChiLog::GetInstance();

  unsigned int last_unknown_end = -1;
  if (not unknowns.empty())
    last_unknown_end = unknowns.back()->GetMapEnd();

  unsigned int new_unknown_index = unknowns.size();

  if (unk_type == UnknownType::SCALAR)
  {
    unknowns.push_back(new ScalarUnknown(last_unknown_end+1));
  }
  else if (unk_type == UnknownType::VECTOR_2)
  {
    unknowns.push_back(new Vector2Unknown(last_unknown_end+1));
  }
  else if (unk_type == UnknownType::VECTOR_3)
  {
    unknowns.push_back(new Vector3Unknown(last_unknown_end+1));
  }
  else if (unk_type == UnknownType::VECTOR_N)
  {
    if (dimension == 0)
    {
      log.Log(LOG_ALLERROR)
        << "UnknownManager: When adding unknown of type VECTOR_N, "
        << "the dimension must not be 0.";
      exit(EXIT_FAILURE);
    }

    unknowns.push_back(new VectorNUnknown(last_unknown_end+1,dimension));
  }
  else if (unk_type == UnknownType::TENSOR)
  {
    if (dimension == 0 or dimension == 1)
    {
      log.Log(LOG_ALLERROR)
        << "UnknownManager: When adding unknown of type TENSOR, "
        << "the dimension must not be 0 or 1.";
      exit(EXIT_FAILURE);
    }

    //unknowns.push_back(new VectorNUnknown(last_unknown_end+1,dimension));
  }
  else
  {
    log.Log(LOG_ALLERROR)
      << "UnknownManager: Invalid call to AddUnknown(). "
      << "Unknown type is probably not supported yet.";
    exit(EXIT_FAILURE);
  }

  return new_unknown_index;
}

//###################################################################
/**Maps the unknown's component within the storage of a node.*/
unsigned int chi_math::UnknownManager::
  MapUnknown(unsigned int unknown_id, unsigned int component)
{
  auto& log = ChiLog::GetInstance();

  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    log.Log(LOG_ALLERROR)
      << "UnknownManager failed call to MapUnknown";
    exit(EXIT_FAILURE);
  }
  return unknowns[unknown_id]->GetMap(component);
}

//###################################################################
/**Determines the total number of components over all unknowns.*/
unsigned int chi_math::UnknownManager::GetTotalUnknownSize()
{
  if (unknowns.empty())
    return 0;

  return unknowns.back()->GetMapEnd()+1;
}