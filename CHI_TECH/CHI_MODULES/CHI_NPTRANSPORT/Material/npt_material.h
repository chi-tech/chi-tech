#ifndef _npt_material_h
#define _npt_material_h

#include <iostream>
#include <vector>



/** Stores the group to group transfer per moment*/
struct NPT_MATERIAL_MOMENT_TRANFER
{
  std::vector<double*> sigma_sm_gp_to_g;
};

/** Stores the material component data.*/
struct NPT_MATERIAL_COMPONENT
{
  std::string name;
  std::vector<double>  group_boundaries;
  std::vector<double>  sigma_t;
  std::vector<double>  sigma_s;
  std::vector<double>  sigma_a;
  std::vector<NPT_MATERIAL_MOMENT_TRANFER*> sigma_sm_matrices;
};

//###################################################################
/** Contains material components and basic operations.*/
class NPT_MATERIAL
{
public:
  int         id;
  std::string name;
  std::vector<NPT_MATERIAL_COMPONENT*> components;
  std::vector<double>  group_boundaries;
  std::vector<double>  sigma_t;
  std::vector<double>  sigma_s;
  std::vector<double>  sigma_a;
  std::vector<NPT_MATERIAL_MOMENT_TRANFER*> sigma_sm_matrices;

public:
  void SetAsPureAbsorber(int number_of_groups,int L,double sigma_t_in);
};


#endif