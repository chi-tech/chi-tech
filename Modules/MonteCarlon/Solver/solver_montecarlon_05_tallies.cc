#include "solver_montecarlon.h"

#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>

#include <FiniteVolume/CellViews/fv_slab.h>
#include <FiniteVolume/CellViews/fv_polygon.h>

typedef unsigned long long TULL;

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog chi_log;
extern ChiMPI chi_mpi;

//###################################################################
/**Makes a contribution to tallies*/
void chi_montecarlon::Solver::ContributeTally(chi_montecarlon::Particle *prtcl,
                                              chi_mesh::Vector pf)
{
  auto cell = grid->cells[prtcl->cur_cell_ind];
  int cell_local_ind = cell->cell_local_id;

  int ir = cell_local_ind*num_grps + prtcl->egrp;

  double tracklength = (pf - prtcl->pos).Norm();

  phi_tally_contrib[ir] += (tracklength*prtcl->w);

  if (std::isnan(tracklength))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Tracklength corruption."
      << " pos  " << prtcl->pos.PrintS()
      << " posf " << pf.PrintS();
    exit(EXIT_FAILURE);
  }

}


//###################################################################
/**Compute tally square contributions.*/
void chi_montecarlon::Solver::ComputeTallySqr()
{
  int num_cells = grid->local_cell_glob_indices.size();
  for (int lc=0; lc<num_cells; lc++)
  {
    int cell_glob_index = grid->local_cell_glob_indices[lc];
    auto cell = grid->cells[cell_glob_index];

    double V = 1.0;
    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      auto cell_fv_view =
        (SlabFVView*)fv_discretization->MapFeView(cell->cell_global_id);
      V = cell_fv_view->volume;
    }
    else if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      auto cell_fv_view =
        (PolygonFVView*)fv_discretization->MapFeView(cell->cell_global_id);
      V = cell_fv_view->volume;
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Unsupported cell type encountered in call to "
        << "chi_montecarlon::Solver::ComputeTallySqr.";
      exit(EXIT_FAILURE);
    }

    int hi = 0;
    int lo = num_grps-1;
    if (group_hi_bound >= 0)
      hi = group_hi_bound;
    if (group_lo_bound >=0 && group_lo_bound >= group_hi_bound)
      lo = group_lo_bound;

    for (int g=hi; g<=lo; g++)
    {
      int ir = lc*num_grps + g;

      phi_tally[ir]     += phi_tally_contrib[ir]/V;
      phi_tally_sqr[ir] += phi_tally_contrib[ir]*phi_tally_contrib[ir]/V/V;
      phi_tally_contrib[ir] = 0.0;
    }
  }
}



//###################################################################
void chi_montecarlon::Solver::ComputeRelativeStdDev()
{
  max_relative_error = 0.0;

  int num_cells = grid->local_cell_glob_indices.size();
  for (int lc=0; lc<num_cells; lc++)
  {
    //=======
    // PROCESS CELL VOLUME HERE
    //=======



    int hi = 0;
    int lo = num_grps-1;
    if (group_hi_bound >= 0)
      hi = group_hi_bound;
    if (group_lo_bound >=0 && group_lo_bound >= group_hi_bound)
      lo = group_lo_bound;

    for (int g=0; g<num_grps; g++)
    {
      int ir = lc*num_grps + g;

      TULL n = batch_sizes[current_batch];

      double x_avg  = phi_global[ir]/nps_global;
      double x2_avg = (phi_global_tally_sqr[ir]/n);

      double stddev = sqrt((x2_avg - x_avg*x_avg)/nps_global);

      if (!std::isinf(stddev/x_avg))
      {
        phi_local_relsigma[ir] = stddev/x_avg;

        if (g==0 && lc == 0)
        {
          max_relative_error = phi_local_relsigma[ir];
          max_relative_error2 = x2_avg;
          max_relative_error3 = x_avg*x_avg;
        }
      }
      else
      {
        phi_local_relsigma[ir] = 0.0;
      }


      if (std::isinf(phi_local_relsigma[ir]))
      {
        printf("Infinite lc=%d g=%d\n", lc,g);
        chi_log.Log(LOG_ALL) << "stddev=" << stddev << "\n";
        chi_log.Log(LOG_ALL) << "xavg=" << x_avg << "\n";
        chi_log.Log(LOG_ALL) << "x2avg=" << x2_avg << "\n";
        chi_log.Log(LOG_ALL) << "phirel=" << phi_local_relsigma[ir] << "\n";
        chi_log.Log(LOG_ALL) << "batch=" << n << "\n";
        chi_log.Log(LOG_ALL) << "nps=" << nps_global << "\n";
      }
    }//for g
  }//for local cell
}


//###################################################################
/**Merges tallies from multiple locations.*/
void chi_montecarlon::Solver::RendesvouzTallies()
{
  TULL temp_nps_global = 0;
  MPI_Allreduce(&nps,&temp_nps_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);

  nps_global += temp_nps_global;

  //============================================= Merge phi_local
  int num_values = phi_tally.size();
  std::vector<double> temp_phi_global(num_values,0.0);

  MPI_Allreduce(phi_tally.data(),
                temp_phi_global.data(),
                num_values,MPI_DOUBLE,
                MPI_SUM,MPI_COMM_WORLD);

  for (int i=0; i<num_values; i++)
    phi_global[i] += temp_phi_global[i];

  //============================================= Merge phi_tally local
  phi_global_tally_sqr.assign(num_values,0.0);

  MPI_Allreduce(phi_tally_sqr.data(),
                phi_global_tally_sqr.data(),
                num_values,MPI_DOUBLE,
                MPI_SUM,MPI_COMM_WORLD);

  //============================================= Reset tallies
  phi_tally.assign(phi_tally.size(),0.0);
  phi_tally_sqr.assign(phi_tally_sqr.size(),0.0);
}
