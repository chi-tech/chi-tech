#include "lbs_linear_boltzmann_solver.h"

#include "FiniteElement/PiecewiseLinear/pwl.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include <iomanip>

//###################################################################
/**Compute balance.*/
void LinearBoltzmann::Solver::ComputeBalance()
{
  auto pwld =
    std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(discretization);
  if (not pwld) throw std::logic_error("Trouble getting PWLD-SDM in " +
                                      std::string(__FUNCTION__));
  SpatialDiscretization_PWLD& grid_fe_view = *pwld;

  //======================================== Zero outflow
  for (auto& cell_transport_view : cell_transport_views)
    cell_transport_view.ZeroOutflow();

  //======================================== Sweep all groupsets to populate
  //                                         outflow
  auto phi_temp = phi_new_local;
  for (auto& groupset : group_sets)
  {
    auto sweep_chunk = SetSweepChunk(groupset);
    MainSweepScheduler sweep_scheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
                                       groupset.angle_agg,
                                       *sweep_chunk);


    sweep_scheduler.sweep_chunk.SetDestinationPhi(phi_temp);
    sweep_scheduler.sweep_chunk.SetSurfaceSourceActiveFlag(true);
    SetSource(groupset,APPLY_MATERIAL_SOURCE);

    phi_temp.assign(phi_temp.size(),0.0);
    sweep_scheduler.Sweep();
  }//for groupset

  //======================================== Compute absorbtion, material-source
  //                                         and in-flow
  double out_flow=0.0;
  double in_flow=0.0;
  double absorbtion=0.0;
  double mat_source=0.0;
  size_t num_groups=groups.size();
  for (auto& cell : grid->local_cells)
  {
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& fe_intgrl_values = grid_fe_view.GetUnitIntegrals(cell);
    const size_t num_nodes = transport_view.NumNodes();

    const auto& IntV_shapeI = fe_intgrl_values.GetIntV_shapeI();

    for (int g=0; g<num_groups; ++g)
      out_flow += transport_view.GetOutflow(g);

    auto& sigma_ag = material_xs[transport_view.XSMapping()]->sigma_ag;

    for (int i=0; i<num_nodes; ++i)
      for (int g=0; g<num_groups; ++g)
      {
        size_t imap = transport_view.MapDOF(i,0,g);
        double phi_0g = phi_old_local[imap];
        double q_0g   = q_moments_local[imap];

        absorbtion += sigma_ag[g]*phi_0g*IntV_shapeI[i];
        mat_source += q_0g*IntV_shapeI[i];
      }//for g
  }//for cell

  //======================================== Consolidate local balances
  double local_balance = mat_source + in_flow - absorbtion - out_flow;
  double globl_balance = 0.0;

  MPI_Allreduce(&local_balance,   //sendbuf
                &globl_balance,   //recvbuf
                1,MPI_DOUBLE,     //count + datatype
                MPI_SUM,          //operation
                MPI_COMM_WORLD);  //communicator

  chi_log.Log() << "Global balance: "
                << std::setprecision(6) << std::scientific
                << globl_balance;
}