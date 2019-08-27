#include "lbs_linear_boltzman_solver.h"
#include <CHI_MESH/CHI_CELL/cell.h>
#include <CHI_MESH/CHI_CELL/cell_slab.h>
#include <CHI_MESH/CHI_CELL/cell_polygon.h>
#include <CHI_MESH/CHI_CELL/cell_polyhedron.h>


//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void CHI_NPTRANSPORT::
AssembleVector(NPT_GROUPSET *groupset, Vec x, double *y)
{
  double* x_ref;
  VecGetArray(x,&x_ref);

  int gsi = groupset->groups[0]->id;
  int gsf = groupset->groups.back()->id;
  int gss = gsf-gsi+1;

  int index = -1;
  for (int c=0; c<grid->local_cell_glob_indices.size(); c++)
  {
    int cell_g_index = grid->local_cell_glob_indices[c];
    auto cell        = grid->cells[cell_g_index];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      for (int i=0; i<2; i++)
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
    }//if slab

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      chi_mesh::CellPolygon* poly_cell =
        (chi_mesh::CellPolygon*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      for (int i=0; i<poly_cell->v_indices.size(); i++)
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
    }//if polygon

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      chi_mesh::CellPolyhedron* polyh_cell =
        (chi_mesh::CellPolyhedron*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      for (int i=0; i<polyh_cell->v_indices.size(); i++)
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
    }//if polyhedron

  }//for cell

  VecRestoreArray(x,&x_ref);
}

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void CHI_NPTRANSPORT::
DisAssembleVector(NPT_GROUPSET *groupset, Vec x_src, double *y)
{
  const double* x_ref;
  VecGetArrayRead(x_src,&x_ref);

  int gsi = groupset->groups[0]->id;
  int gsf = groupset->groups.back()->id;
  int gss = gsf-gsi+1;

  int index = -1;
  for (int c=0; c<grid->local_cell_glob_indices.size(); c++)
  {
    int cell_g_index = grid->local_cell_glob_indices[c];
    auto cell        = grid->cells[cell_g_index];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      for (int i=0; i<2; i++)
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
    }//if slab

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      chi_mesh::CellPolygon* poly_cell =
        (chi_mesh::CellPolygon*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      for (int i=0; i<poly_cell->v_indices.size(); i++)
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
    }//if polygon

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      chi_mesh::CellPolyhedron* polyh_cell =
        (chi_mesh::CellPolyhedron*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      for (int i=0; i<polyh_cell->v_indices.size(); i++)
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
    }//if polyhedron

  }//for cell

  VecRestoreArrayRead(x_src,&x_ref);
}


//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void CHI_NPTRANSPORT::
DisAssembleVectorLocalToLocal(NPT_GROUPSET *groupset, double* x_src, double *y)
{
  const double* x_ref=x_src;

  int gsi = groupset->groups[0]->id;
  int gsf = groupset->groups.back()->id;
  int gss = groupset->groups.size();

  int index = -1;
  for (int c=0; c<grid->local_cell_glob_indices.size(); c++)
  {
    int cell_g_index = grid->local_cell_glob_indices[c];
    auto cell        = grid->cells[cell_g_index];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      for (int i=0; i<2; i++)
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
    }//if slab

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      chi_mesh::CellPolygon* poly_cell =
        (chi_mesh::CellPolygon*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      for (int i=0; i<poly_cell->v_indices.size(); i++)
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
    }//if polygon

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      chi_mesh::CellPolyhedron* polyh_cell =
        (chi_mesh::CellPolyhedron*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      for (int i=0; i<polyh_cell->v_indices.size(); i++)
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
    }//if polyhedron


  }//for cell

}