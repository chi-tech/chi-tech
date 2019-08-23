#ifndef _chi_npt_structs_h
#define _chi_npt_structs_h

#include "../../../CHI_RESOURCES/Dependencies/Eigen/Dense"


#define PARTITION_METHOD_SERIAL        1
#define PARTITION_METHOD_FROM_SURFACE  2

/**Struct for storing NPT options.*/
struct NPT_OPTIONS
{
  int scattering_order;
  int partition_method;
  int sweep_eager_limit;

  NPT_OPTIONS()
  {
    scattering_order = 0;
    partition_method = PARTITION_METHOD_SERIAL;
    sweep_eager_limit= 32000;
  }
};

class NPT_CELLVIEW
{
public:
  std::vector<int> node_dof_mapping;
};



/**Transport view of a cell*/
class NPT_CELLVIEW_FULL : public NPT_CELLVIEW
{
public:
  std::vector<bool>  face_f_upwind_flag;
  std::vector<int>   face_f_adj_part_id;
  std::vector<int>   face_boundary_id;

  int dof_phi_map_start;
  int dofs;
  int xs_id;

private:
  int num_grps;
  int num_moms;

public:
  NPT_CELLVIEW_FULL(int in_dofs,int num_G, int num_m)
  {
    dof_phi_map_start = -1;
    dofs = in_dofs;
    num_grps = num_G;
    num_moms = num_m;
  }

  int MapDOF(int dof, int moment, int grp)
  {
    return dof_phi_map_start + dof*num_grps*num_moms + num_grps*moment + grp;
  }
};



#endif