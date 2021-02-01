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
    last_unknown_end = unknowns.back().GetMapEnd();

  unsigned int new_unknown_index = unknowns.size();

  if (unk_type == UnknownType::SCALAR)
  {
    unknowns.emplace_back(UnknownType::SCALAR,1,last_unknown_end+1);
    unknowns.back().text_name = "Unknown_" + std::to_string(unknowns.size()-1);
  }
  else if (unk_type == UnknownType::VECTOR_2)
  {
    unknowns.emplace_back(UnknownType::VECTOR_2,2,last_unknown_end+1);
    unknowns.back().text_name = "Unknown_" + std::to_string(unknowns.size()-1);
  }
  else if (unk_type == UnknownType::VECTOR_3)
  {
    unknowns.emplace_back(UnknownType::VECTOR_3,3,last_unknown_end+1);
    unknowns.back().text_name = "Unknown_" + std::to_string(unknowns.size()-1);
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

    unknowns.emplace_back(UnknownType::VECTOR_N,dimension,last_unknown_end+1);
    unknowns.back().text_name = "Unknown_" + std::to_string(unknowns.size()-1);
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

    throw std::invalid_argument("UnknownManager: TENSOR unknowns are not "
                                "supported yet.");
  }
  else
  {
    throw std::logic_error("UnknownManager: Invalid call to AddUnknown(). "
                           "Unknown type is probably not supported yet.");
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
  return unknowns[unknown_id].GetMap(component);
}

//###################################################################
/**Determines the total number of components over all unknowns.*/
unsigned int chi_math::UnknownManager::GetTotalUnknownSize()
{
  if (unknowns.empty())
    return 0;

  return unknowns.back().GetMapEnd()+1;
}

//###################################################################
/**Sets the number of off block connections for the given unknown.
 * All the components will be set to the same amount.*/
void chi_math::UnknownManager::
  SetUnknownNumOffBlockConnections(unsigned int unknown_id,
                                   int num_conn)
{
  auto& log = ChiLog::GetInstance();

  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    log.Log(LOG_ALLERROR)
      << "UnknownManager failed call to SetUnknownNumOffBlockConnections,"
         " illegal index. " << unknown_id;
    exit(EXIT_FAILURE);
  }

  for (auto& val : unknowns[unknown_id].num_off_block_connections)
    val = num_conn;
}
//###################################################################
/**Sets the number of off block connections for the given unknown-
 * component pair.*/
void chi_math::UnknownManager::
  SetComponentNumOffBlockConnections(unsigned int unknown_id,
                                     unsigned int component,
                                     int num_conn)
{
  auto& log = ChiLog::GetInstance();

  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    log.Log(LOG_ALLERROR)
      << "UnknownManager failed call to SetUnknownComponentTextName,"
         " illegal unknown index. " << unknown_id;
    exit(EXIT_FAILURE);
  }

  if (component < 0 or component >= unknowns[unknown_id].num_components)
  {
    log.Log(LOG_ALLERROR)
      << "UnknownManager failed call to SetUnknownComponentTextName,"
         " illegal component index. " << component;
    exit(EXIT_FAILURE);
  }

  unknowns[unknown_id].num_off_block_connections[component] = num_conn;

}

//###################################################################
/**Sets a text name for the indicated unknown.*/
void chi_math::UnknownManager::
  SetUnknownTextName(unsigned int unknown_id,
                     const std::string& in_text_name)
{
  auto& log = ChiLog::GetInstance();

  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    log.Log(LOG_ALLERROR)
      << "UnknownManager failed call to SetUnknownTextName,"
         " illegal index. " << unknown_id;
    exit(EXIT_FAILURE);
  }

  unknowns[unknown_id].text_name = in_text_name;
}

//###################################################################
/**Sets the text name to be associated with each component of the
 * unknown.*/
void chi_math::UnknownManager::
  SetUnknownComponentTextName(unsigned int unknown_id,
                              unsigned int component,
                              const std::string& in_text_name)
{
  auto& log = ChiLog::GetInstance();

  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    log.Log(LOG_ALLERROR)
      << "UnknownManager failed call to SetUnknownComponentTextName,"
         " illegal unknown index. " << unknown_id;
    exit(EXIT_FAILURE);
  }

  if (component < 0 or component >= unknowns[unknown_id].num_components)
  {
    log.Log(LOG_ALLERROR)
      << "UnknownManager failed call to SetUnknownComponentTextName,"
         " illegal component index. " << component;
    exit(EXIT_FAILURE);
  }

  unknowns[unknown_id].component_text_names[component] = in_text_name;

}