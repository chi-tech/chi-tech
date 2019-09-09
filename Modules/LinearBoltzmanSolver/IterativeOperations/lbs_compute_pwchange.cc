#include "../lbs_linear_boltzman_solver.h"
#include <ChiMesh/Cell/cell.h>
#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>

//###################################################################
/**Computes the point wise change between phi_new and phi_old.*/
double LinearBoltzmanSolver::ComputePiecewiseChange(LBSGroupset* groupset)
{
  double pw_change = 0.0;
  double sum_m0 = 0.0;

  int gsi = groupset->groups[0]->id;
  int gsf = groupset->groups.back()->id;
  int deltag = groupset->groups.size();

  for (int c=0; c<grid->local_cell_glob_indices.size(); c++)
  {
    int cell_g_index = grid->local_cell_glob_indices[c];
    auto cell        = grid->cells[cell_g_index];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      chi_mesh::CellSlab* slab_cell =
        (chi_mesh::CellSlab*)cell;
      LBSCellViewFull* transport_view =
        (LBSCellViewFull*)cell_transport_views[c];

      for (int i=0; i<2; i++)
      {
        for (int m=0; m<num_moments; m++)
        {
          int mapping = transport_view->MapDOF(i,m,gsi);
          double* phi_new_m = &phi_new_local.data()[mapping];
          double* phi_old_m = &phi_old_local.data()[mapping];

          for (int g=0; g<deltag; g++)
          {
            int map0 = transport_view->MapDOF(i,0,gsi+g);

            double abs_phi_m0     = fabs(phi_new_local[map0]);
            double abs_phi_old_m0 = fabs(phi_old_local[map0]);
            double max_phi = std::max(abs_phi_m0,abs_phi_old_m0);

            double delta_phi = std::fabs(phi_new_m[g] - phi_old_m[g]);

            if (max_phi >= std::numeric_limits<double>::min())
              pw_change = std::max(delta_phi/max_phi,pw_change);
            else
              pw_change = std::max(delta_phi,pw_change);

          }//for g
        }//for m
      }//for i
    }//if slab

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      chi_mesh::CellPolygon* poly_cell =
        (chi_mesh::CellPolygon*)cell;
      LBSCellViewFull* transport_view =
        (LBSCellViewFull*)cell_transport_views[c];

      for (int i=0; i<poly_cell->v_indices.size(); i++)
      {
        for (int m=0; m<num_moments; m++)
        {
          int mapping = transport_view->MapDOF(i,m,gsi);
          double* phi_new_m = &phi_new_local.data()[mapping];
          double* phi_old_m = &phi_old_local.data()[mapping];

          for (int g=0; g<deltag; g++)
          {
            int map0 = transport_view->MapDOF(i,0,gsi+g);

            double abs_phi_m0     = fabs(phi_new_local[map0]);
            double abs_phi_old_m0 = fabs(phi_old_local[map0]);
            double max_phi = std::max(abs_phi_m0,abs_phi_old_m0);

            double delta_phi = std::fabs(phi_new_m[g] - phi_old_m[g]);

            if (max_phi >= std::numeric_limits<double>::min())
              pw_change = std::max(delta_phi/max_phi,pw_change);
            else
              pw_change = std::max(delta_phi,pw_change);

          }//for g
        }//for m
      }//for i
    }//if polyhedron

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      chi_mesh::CellPolyhedron* polyh_cell =
        (chi_mesh::CellPolyhedron*)cell;
      LBSCellViewFull* transport_view =
        (LBSCellViewFull*)cell_transport_views[c];

      for (int i=0; i<polyh_cell->v_indices.size(); i++)
      {
        for (int m=0; m<num_moments; m++)
        {
          int mapping = transport_view->MapDOF(i,m,gsi);
          double* phi_new_m = &phi_new_local.data()[mapping];
          double* phi_old_m = &phi_old_local.data()[mapping];

          for (int g=0; g<deltag; g++)
          {
            int map0 = transport_view->MapDOF(i,0,gsi+g);

            double abs_phi_m0     = fabs(phi_new_local[map0]);
            double abs_phi_old_m0 = fabs(phi_old_local[map0]);
            double max_phi = std::max(abs_phi_m0,abs_phi_old_m0);

            double delta_phi = std::fabs(phi_new_m[g] - phi_old_m[g]);

            if (max_phi >= std::numeric_limits<double>::min())
              pw_change = std::max(delta_phi/max_phi,pw_change);
            else
              pw_change = std::max(delta_phi,pw_change);

          }//for g
        }//for m
      }//for i
    }//if polyhedron

  }//for c

//  const real8 abs_phi_0 = fabs(phi_0);
//  const real8 abs_old_phi_0 = fabs(old_phi_0);
//  const real8 maxv = std::max(abs_phi_0, abs_old_phi_0);
//  const real8 diff = fabs(phi - old_phi);
//  const real8 pw_change = (maxv >= std::numeric_limits<real8>::min()) ?
//                          (diff / maxv) : diff;

  double global_pw_change = 0.0;

  MPI_Allreduce(&pw_change,&global_pw_change,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  return global_pw_change;
}