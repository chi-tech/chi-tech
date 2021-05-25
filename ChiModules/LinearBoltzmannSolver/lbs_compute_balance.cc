#include "lbs_linear_boltzmann_solver.h"

#include "FiniteElement/PiecewiseLinear/pwl.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include <iomanip>

//###################################################################
/**Compute balance.*/
void LinearBoltzmann::Solver::ComputeBalance()
{
  chi_log.Log() << "Computing balance.";

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
  auto mat_src = phi_temp;
  mat_src.assign(mat_src.size(),0.0);
  int gs=0;
  for (auto& groupset : group_sets)
  {
    chi_log.Log() << "******************** Sweeping GS " << gs;
    ComputeSweepOrderings(groupset);
    InitFluxDataStructures(groupset);

    auto sweep_chunk = SetSweepChunk(groupset);
    MainSweepScheduler sweep_scheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
                                       groupset.angle_agg,
                                       *sweep_chunk);

    sweep_scheduler.sweep_chunk.SetDestinationPhi(phi_temp);
    sweep_scheduler.sweep_chunk.SetSurfaceSourceActiveFlag(true);
    SetSource(groupset,APPLY_MATERIAL_SOURCE);

    for (size_t i=0; i<mat_src.size(); ++i)
      mat_src[i] += q_moments_local[i];

    phi_temp.assign(phi_temp.size(),0.0);
    sweep_scheduler.Sweep();

    ResetSweepOrderings(groupset);

    MPI_Barrier(MPI_COMM_WORLD);
    ++gs;
  }//for groupset

  chi_log.Log() << "Computing items";

  //======================================== Compute absorbtion, material-source
  //                                         and in-flow
  double out_flow=0.0;
  double in_flow=0.0;
  double absorbtion=0.0;
  double IntV_q=0.0;
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
        double q_0g   = mat_src[imap];

//        absorbtion += sigma_ag[g]*phi_0g*IntV_shapeI[i];
        IntV_q += q_0g * IntV_shapeI[i];
      }//for g
  }//for cell

  //======================================== Consolidate local balances
  double local_balance = IntV_q + in_flow - absorbtion - out_flow;
  double globl_balance = 0.0;

  MPI_Allreduce(&local_balance,   //sendbuf
                &globl_balance,   //recvbuf
                1,MPI_DOUBLE,     //count + datatype
                MPI_SUM,          //operation
                MPI_COMM_WORLD);  //communicator

  chi_log.Log(LOG_ALL) << "Local balance: "
                << std::setprecision(6) << std::scientific
                << local_balance;
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log() << "Global balance: "
                << std::setprecision(6) << std::scientific
                << globl_balance;
}