#include "lbs_linear_boltzmann_solver.h"
#include <ChiMesh/Cell/cell.h>

//###################################################################
/**Maps the local storage location of phi given a cell, the node of
 * the cell, the moment and the group.*/
int LinearBoltzmann::Solver::
MapDOF(chi_mesh::Cell *cell, int dof, int mom, int g)
{
  int G = groups.size();
  int address = local_cell_dof_array_address[cell->local_id];
//  int block_address = (address+dof)*num_moments*groups.size();

  return (address+dof)*num_moments*G + G*mom + g;
}

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void LinearBoltzmann::Solver::
AssembleVector(LBSGroupset *groupset, Vec x, double *y)
{
  double* x_ref;
  VecGetArray(x,&x_ref);

  int gsi = groupset->groups[0]->id;
  int gsf = groupset->groups.back()->id;
  int gss = gsf-gsi+1;

  int index = -1;
  for (const auto& cell : grid->local_cells)
  {
    auto transport_view =
      (LinearBoltzmann::CellViewFull*)cell_transport_views[cell.local_id];

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      for (int m=0; m<num_moments; m++)
      {
        int mapping = transport_view->MapDOF(i,m,gsi);
        double* source_mapped = &y[mapping];
        for (int g=0; g<gss; g++)
        {
          index++;
          x_ref[index] = source_mapped[g]; //On purpose
        }//for g
      }//for moment
    }//for dof
  }//for cell


  groupset->angle_agg->AssembleAngularUnknowns(index,x_ref);

  VecRestoreArray(x,&x_ref);
}

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void LinearBoltzmann::Solver::
DisAssembleVector(LBSGroupset *groupset, Vec x_src, double *y)
{
  const double* x_ref;
  VecGetArrayRead(x_src,&x_ref);

  int gsi = groupset->groups[0]->id;
  int gsf = groupset->groups.back()->id;
  int gss = gsf-gsi+1;

  int index = -1;
  for (const auto& cell : grid->local_cells)
  {
    auto transport_view =
      (LinearBoltzmann::CellViewFull*)cell_transport_views[cell.local_id];

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      for (int m=0; m<num_moments; m++)
      {
        int mapping = transport_view->MapDOF(i,m,gsi);
        double* destination_mapped = &y[mapping];
        for (int g=0; g<gss; g++)
        {
          index++;
          destination_mapped[g] = x_ref[index];
        }//for g
      }//for moment
    }//for dof
  }//for cell

  groupset->angle_agg->DisassembleAngularUnknowns(index,x_ref);

  VecRestoreArrayRead(x_src,&x_ref);
}


//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void LinearBoltzmann::Solver::
DisAssembleVectorLocalToLocal(LBSGroupset *groupset, double* x_src, double *y)
{
  const double* x_ref=x_src;

  int gsi = groupset->groups[0]->id;
  int gsf = groupset->groups.back()->id;
  int gss = groupset->groups.size();

  int index = -1;
  for (const auto& cell : grid->local_cells)
  {
    auto transport_view =
      (LinearBoltzmann::CellViewFull*)cell_transport_views[cell.local_id];

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      for (int m=0; m<num_moments; m++)
      {
        int mapping = transport_view->MapDOF(i,m,gsi);
        double* destination_mapped = &y[mapping];
        const double* source_mapped      = &x_ref[mapping];
        for (int g=0; g<gss; g++)
        {
          index++;
          destination_mapped[g] = source_mapped[g];
        }//for g
      }//for moment
    }//for dof
  }//for cell

}