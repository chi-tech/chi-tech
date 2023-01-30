#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "LBSSteadyState/Groupset/lbs_groupset.h"

;

#include <iomanip>

//###################################################################
/**Zeroes all the outflow data-structures required to compute
 * balance.*/
void lbs::SteadyStateSolver::ZeroOutflowBalanceVars(LBSGroupset& groupset)
{
  for (auto& cell_transport_view : cell_transport_views)
    for (auto& group : groupset.groups)
      cell_transport_view.ZeroOutflow(group.id);
}

//###################################################################
/**Compute balance.*/
void lbs::SteadyStateSolver::ComputeBalance()
{
  MPI_Barrier(MPI_COMM_WORLD);
  chi::log.Log() << "\n********** Computing balance\n";

  auto pwld =
    std::dynamic_pointer_cast<chi_math::SpatialDiscretization_PWLD>(discretization);
  if (not pwld) throw std::logic_error("Trouble getting PWLD-SDM in " +
                                      std::string(__FUNCTION__));
  chi_math::SpatialDiscretization_PWLD& grid_fe_view = *pwld;

  //======================================== Get material source
  // This is done using the SetSource routine
  // because it allows a lot of flexibility.
  auto mat_src = phi_old_local;
  mat_src.assign(mat_src.size(),0.0);
  for (auto& groupset : groupsets)
  {
    q_moments_local.assign(q_moments_local.size(), 0.0);
    SetSource(groupset, q_moments_local,
              APPLY_FIXED_SOURCES | APPLY_AGS_FISSION_SOURCES |
              APPLY_WGS_FISSION_SOURCES);
    ScopedCopySTLvectors(groupset, q_moments_local, mat_src);
  }

  //======================================== Initialize diffusion params
  //                                         for xs
  // This populates sigma_a
  for (const auto& mat_id_xs : matid_to_xs_map)
  {
    const auto& xs = mat_id_xs.second;
    if (not xs->diffusion_initialized)
      xs->ComputeDiffusionParameters();
  }

  //======================================== Compute absorption, material-source
  //                                         and in-flow
  double local_out_flow   = 0.0;
  double local_in_flow    = 0.0;
  double local_absorption = 0.0;
  double local_production = 0.0;
  for (const auto& cell : grid->local_cells)
  {
    const auto&  transport_view   = cell_transport_views[cell.local_id];
    const auto&  fe_intgrl_values = grid_fe_view.GetUnitIntegrals(cell);
    const size_t num_nodes        = transport_view.NumNodes();
    const auto&  IntV_shapeI      = fe_intgrl_values.GetIntV_shapeI();
    const auto&  IntS_shapeI      = fe_intgrl_values.GetIntS_shapeI();

    //====================================== Inflow
    // This is essentially an integration over
    // all faces, all angles, and all groups.
    // Only the cosines that are negative are
    // added to the integral.
    for (int f=0; f<cell.faces.size(); ++f)
    {
      const auto& face  = cell.faces[f];

      for (const auto& groupset : groupsets)
      {
        for (int n = 0; n < groupset.quadrature->omegas.size(); ++n)
        {
          const auto &omega = groupset.quadrature->omegas[n];
          const double wt = groupset.quadrature->weights[n];
          const double mu = omega.Dot(face.normal);

          if (mu < 0.0 and (not face.has_neighbor)) //mu<0 and bndry
          {
            const auto &bndry = sweep_boundaries[face.neighbor_id];
            for (int fi = 0; fi < face.vertex_ids.size(); ++fi)
            {
              const int i = fe_intgrl_values.FaceDofMapping(f, fi);
              const auto &IntFi_shapeI = IntS_shapeI[f][i];

              for (const auto &group : groupset.groups)
              {
                const int g = group.id;
                const double psi = bndry->boundary_flux[g];
                local_in_flow -= mu * wt * psi * IntFi_shapeI;
              }//for g
            }//for fi
          }//if bndry
        }//for n
      }//for groupset
    }//for f


    //====================================== Outflow
    //The group-wise outflow was determined
    //during a solve so here we just
    //consolidate it.
    for (int g=0; g<num_groups; ++g)
      local_out_flow += transport_view.GetOutflow(g);

    //====================================== Absorption and Src
    //Isotropic flux based absorption and source
    auto& xs = transport_view.XS();
    for (int i=0; i<num_nodes; ++i)
      for (int g=0; g<num_groups; ++g)
      {
        size_t imap   = transport_view.MapDOF(i,0,g);
        double phi_0g = phi_old_local[imap];
        double q_0g   = mat_src[imap];

        local_absorption += xs.sigma_a[g] * phi_0g * IntV_shapeI[i];
        local_production += q_0g * IntV_shapeI[i];
      }//for g
  }//for cell

  //======================================== Consolidate local balances
  double local_balance = local_production + local_in_flow
                       - local_absorption - local_out_flow;
  double local_gain    = local_production + local_in_flow;

  std::vector<double> local_balance_table = {local_absorption,
                                             local_production,
                                             local_in_flow,
                                             local_out_flow,
                                             local_balance,
                                             local_gain};
  size_t table_size = local_balance_table.size();

  std::vector<double> globl_balance_table(table_size,0.0);

  MPI_Allreduce(local_balance_table.data(),      //sendbuf
                globl_balance_table.data(),      //recvbuf
                table_size,MPI_DOUBLE,           //count + datatype
                MPI_SUM,                         //operation
                MPI_COMM_WORLD);                 //communicator

  double globl_absorption = globl_balance_table.at(0);
  double globl_production = globl_balance_table.at(1);
  double globl_in_flow    = globl_balance_table.at(2);
  double globl_out_flow   = globl_balance_table.at(3);
  double globl_balance    = globl_balance_table.at(4);
  double globl_gain       = globl_balance_table.at(5);

  chi::log.Log() << "Balance table:\n"
    << std::setprecision(5) << std::scientific
    << " Absorption rate          = " << globl_absorption               << "\n"
    << " Production rate          = " << globl_production               << "\n"
    << " In-flow rate             = " << globl_in_flow                  << "\n"
    << " Out-flow rate            = " << globl_out_flow                 << "\n"
    << " Integrated scalar flux   = " << globl_gain                     << "\n"
    << " Net Gain/Loss            = " << globl_balance                  << "\n"
    << " Net Gain/Loss normalized = " << globl_balance/globl_gain       << "\n";

  chi::log.Log() << "\n********** Done computing balance\n";

  MPI_Barrier(MPI_COMM_WORLD);
}