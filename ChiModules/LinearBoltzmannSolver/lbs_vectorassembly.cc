#include "lbs_linear_boltzmann_solver.h"

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void LinearBoltzmann::Solver::
  SetPETScVecFromSTLvector(LBSGroupset& groupset, Vec x,
                           const std::vector<double>& y,
                           bool with_delayed_psi/*=false*/)
{
  double* x_ref;
  VecGetArray(x,&x_ref);

  int gsi = groupset.groups[0].id;
  int gsf = groupset.groups.back().id;
  int gss = gsf-gsi+1;

  int index = -1;
  for (const auto& cell : grid->local_cells)
  {
    auto& transport_view = cell_transport_views[cell.local_id];

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      for (int m=0; m<num_moments; m++)
      {
        size_t mapping = transport_view.MapDOF(i,m,gsi);
        for (int g=0; g<gss; g++)
        {
          index++;
          x_ref[index] = y[mapping+g]; //Offset on purpose
        }//for g
      }//for moment
    }//for dof
  }//for cell

  if (with_delayed_psi)
    groupset.angle_agg.AppendDelayedAngularDOFsToArray(index, x_ref);

  VecRestoreArray(x,&x_ref);
}

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void LinearBoltzmann::Solver::
  SetSTLvectorFromPETScVec(LBSGroupset& groupset, Vec x_src,
                           std::vector<double>& y,
                           bool with_delayed_psi/*=false*/)
{
  const double* x_ref;
  VecGetArrayRead(x_src,&x_ref);

  int gsi = groupset.groups[0].id;
  int gsf = groupset.groups.back().id;
  int gss = gsf-gsi+1;

  int index = -1;
  for (const auto& cell : grid->local_cells)
  {
    auto& transport_view = cell_transport_views[cell.local_id];

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      for (int m=0; m<num_moments; m++)
      {
        size_t mapping = transport_view.MapDOF(i,m,gsi);
        for (int g=0; g<gss; g++)
        {
          index++;
          y[mapping+g] = x_ref[index];
        }//for g
      }//for moment
    }//for dof
  }//for cell

  if (with_delayed_psi)
    groupset.angle_agg.SetDelayedAngularDOFsFromArray(index, x_ref);

  VecRestoreArrayRead(x_src,&x_ref);
}


//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void LinearBoltzmann::Solver::
  ScopedCopySTLvectors(LBSGroupset& groupset,
                       const std::vector<double>& x_src,
                       std::vector<double>& y)
{
  int gsi = groupset.groups[0].id;
  size_t gss = groupset.groups.size();

  int index = -1;
  for (const auto& cell : grid->local_cells)
  {
    auto& transport_view = cell_transport_views[cell.local_id];

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      for (int m=0; m<num_moments; m++)
      {
        size_t mapping = transport_view.MapDOF(i,m,gsi);
        for (int g=0; g<gss; g++)
        {
          index++;
          y[mapping+g] = x_src[mapping+g];
        }//for g
      }//for moment
    }//for dof
  }//for cell

}