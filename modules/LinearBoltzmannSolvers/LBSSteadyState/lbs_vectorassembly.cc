#include "lbs_linear_boltzmann_solver.h"
#include "LBSSteadyState/Groupset/lbs_groupset.h"

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void lbs::SteadyStateSolver::
  SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset, Vec x,
                                    const std::vector<double>& y,
                                    bool with_delayed_psi/*=false*/)
{
  double* x_ref;
  VecGetArray(x,&x_ref);

  int gsi = groupset.groups[0].id;
  int gsf = groupset.groups.back().id;
  int gss = gsf-gsi+1;

  int index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id];

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      for (int m=0; m < num_moments_; m++)
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
void lbs::SteadyStateSolver::
  SetGSSTLvectorFromPrimarySTLvector(LBSGroupset& groupset,
                                     std::vector<double>& x,
                                     const std::vector<double>& y,
                                     bool with_delayed_psi/*=false*/)
{
  int gsi = groupset.groups[0].id;
  int gsf = groupset.groups.back().id;
  int gss = gsf-gsi+1;

  int index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id];

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      for (int m=0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i,m,gsi);
        for (int g=0; g<gss; g++)
        {
          index++;
          x[index] = y[mapping+g]; //Offset on purpose
        }//for g
      }//for moment
    }//for dof
  }//for cell

  if (with_delayed_psi)
    groupset.angle_agg.AppendDelayedAngularDOFsToArray(index, x.data());
}

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void lbs::SteadyStateSolver::
  SetPrimarySTLvectorFromGSPETScVec(LBSGroupset& groupset, Vec x_src,
                                    std::vector<double>& y,
                                    bool with_delayed_psi/*=false*/)
{
  const double* x_ref;
  VecGetArrayRead(x_src,&x_ref);

  int gsi = groupset.groups[0].id;
  int gsf = groupset.groups.back().id;
  int gss = gsf-gsi+1;

  int index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id];

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      for (int m=0; m < num_moments_; m++)
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
void lbs::SteadyStateSolver::
  SetPrimarySTLvectorFromGSSTLvector(LBSGroupset& groupset,
                                     const std::vector<double>& x_src,
                                     std::vector<double>& y,
                                     bool with_delayed_psi/*=false*/)
{
  int gsi = groupset.groups[0].id;
  int gsf = groupset.groups.back().id;
  int gss = gsf-gsi+1;

  int index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id];

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      for (int m=0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i,m,gsi);
        for (int g=0; g<gss; g++)
        {
          index++;
          y[mapping+g] = x_src[index];
        }//for g
      }//for moment
    }//for dof
  }//for cell

  if (with_delayed_psi)
    groupset.angle_agg.SetDelayedAngularDOFsFromArray(index, x_src.data());
}


//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void lbs::SteadyStateSolver::
  GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                const std::vector<double>& x_src,
                                std::vector<double>& y,
                                bool with_delayed_psi/*=false*/)
{
  int gsi = groupset.groups[0].id;
  size_t gss = groupset.groups.size();

  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id];

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      for (int m=0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i,m,gsi);
        for (int g=0; g<gss; g++)
        {
          y[mapping+g] = x_src[mapping+g];
        }//for g
      }//for moment
    }//for dof
  }//for cell

  if (with_delayed_psi)
    groupset.angle_agg.SetDelayedPsiNew2Old();

}