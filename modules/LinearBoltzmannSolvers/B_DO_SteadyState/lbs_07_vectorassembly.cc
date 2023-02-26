#include "lbs_DO_steady_state.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void lbs::DiscOrdSteadyStateSolver::
  SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset, Vec x,
                                    PhiSTLOption which_phi)
{
  const std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW: y_ptr = &phi_new_local_; break;
    case PhiSTLOption::PHI_OLD: y_ptr = &phi_old_local_; break;
    default:
      throw std::logic_error("SetGSPETScVecFromPrimarySTLvector");
  }

  double* x_ref;
  VecGetArray(x,&x_ref);

  int gsi = groupset.groups_.front().id_;
  int gsf = groupset.groups_.back().id_;
  int gss = gsf-gsi+1;

  int index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i=0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m=0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i,m,gsi);
        for (int g=0; g<gss; g++)
        {
          index++;
          x_ref[index] = (*y_ptr)[mapping+g]; //Offset on purpose
        }//for g
      }//for moment
    }//for dof
  }//for cell

  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      groupset.angle_agg_.AppendNewDelayedAngularDOFsToArray(index, x_ref);
      break;
    case PhiSTLOption::PHI_OLD:
      groupset.angle_agg_.AppendOldDelayedAngularDOFsToArray(index, x_ref);
      break;
  }

  VecRestoreArray(x,&x_ref);
}

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void lbs::DiscOrdSteadyStateSolver::
  SetPrimarySTLvectorFromGSPETScVec(LBSGroupset& groupset, Vec x_src,
                                    PhiSTLOption which_phi)
{
  std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW: y_ptr = &phi_new_local_; break;
    case PhiSTLOption::PHI_OLD: y_ptr = &phi_old_local_; break;
    default:
      throw std::logic_error("SetPrimarySTLvectorFromGSPETScVec");
  }

  const double* x_ref;
  VecGetArrayRead(x_src,&x_ref);

  int gsi = groupset.groups_.front().id_;
  int gsf = groupset.groups_.back().id_;
  int gss = gsf-gsi+1;

  int index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i=0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m=0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i,m,gsi);
        for (int g=0; g<gss; g++)
        {
          index++;
          (*y_ptr)[mapping+g] = x_ref[index];
        }//for g
      }//for moment
    }//for dof
  }//for cell

  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      groupset.angle_agg_.SetNewDelayedAngularDOFsFromArray(index, x_ref);
      break;
    case PhiSTLOption::PHI_OLD:
      groupset.angle_agg_.SetOldDelayedAngularDOFsFromArray(index, x_ref);
  }

  VecRestoreArrayRead(x_src,&x_ref);
}

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void lbs::DiscOrdSteadyStateSolver::
  GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                PhiSTLOption from_which_phi,
                                PhiSTLOption to_which_phi)
{
  std::vector<double>* y_ptr;
  switch (to_which_phi)
  {
    case PhiSTLOption::PHI_NEW: y_ptr = &phi_new_local_; break;
    case PhiSTLOption::PHI_OLD: y_ptr = &phi_old_local_; break;
    default:
      throw std::logic_error("GSScopedCopyPrimarySTLvectors");
  }

  std::vector<double>* x_src_ptr;
  switch (from_which_phi)
  {
    case PhiSTLOption::PHI_NEW: x_src_ptr = &phi_new_local_; break;
    case PhiSTLOption::PHI_OLD: x_src_ptr = &phi_old_local_; break;
    default:
      throw std::logic_error("GSScopedCopyPrimarySTLvectors");
  }

  int gsi = groupset.groups_.front().id_;
  size_t gss = groupset.groups_.size();

  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i=0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m=0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i,m,gsi);
        for (int g=0; g<gss; g++)
        {
          (*y_ptr)[mapping+g] = (*x_src_ptr)[mapping+g];
        }//for g
      }//for moment
    }//for dof
  }//for cell

  if (from_which_phi == PhiSTLOption::PHI_NEW and
      to_which_phi == PhiSTLOption::PHI_OLD)
    groupset.angle_agg_.SetDelayedPsiOld2New();
  if (from_which_phi == PhiSTLOption::PHI_OLD and
      to_which_phi == PhiSTLOption::PHI_NEW)
    groupset.angle_agg_.SetDelayedPsiNew2Old();
}