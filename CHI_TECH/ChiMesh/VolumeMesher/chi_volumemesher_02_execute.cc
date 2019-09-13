#include "chi_volumemesher.h"
#include <iostream>



void chi_mesh::VolumeMesher::Execute()
{
  std::cout << "Empty volume mesher, nothing executed.";
  std::cout << std::endl;
}

//###################################################################
/** Maps a node index to a reordered index if it exists.*/
int chi_mesh::VolumeMesher::MapNode(int iref)
{
  //============================================= Check if reordering has
  //                                              been done
  if (this->node_ordering.size()==0)
  {
    return iref;
  }

  //============================================= If it has been done
  chi_mesh::NodeIndexMap* mapping;

  try{
    mapping = this->node_ordering.at(iref);
  }
  catch(const std::out_of_range& o){
    std::cerr << "ERROR: Cannot reference node mapping.";
    exit(EXIT_FAILURE);
  }
//  mapping = this->node_ordering.at(iref);

  return mapping->mapped_to;
//  return node_ordering2[iref];
}

//###################################################################
/** Maps a node index to a reordered index if it exists.*/
int chi_mesh::VolumeMesher::ReverseMapNode(int i)
{
  //============================================= Check if reordering has
  //                                              been done
  if (this->reverse_node_ordering.size()==0)
  {
    return i;
  }

  //============================================= If it has been done
  return reverse_node_ordering[i];
}
