#include "source_function.h"

#include "A_LBSSolver/lbs_solver.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace lbs
{

//###################################################################
/**Constructor.*/
SourceFunction::SourceFunction(const LBSSolver &lbs_solver) :
  lbs_solver_(lbs_solver)
{}


//###################################################################
/**Sets the source moments for the groups in the current group set.
 *
 * \param groupset The groupset the under consideration.
 * \param destination_q A vector to contribute the source to.
 * \param phi_local The primary STL vector to operate off.
 * \param source_flags Flags for adding specific terms into the
 *        destination vector. Available flags are for applying
 *        the material source, across/within-group scattering,
 *        and across/within-groups fission.
 *
 * */
void SourceFunction::operator()(LBSGroupset &groupset,
                                std::vector<double> &destination_q,
                                const std::vector<double> &phi_local,
                                SourceFlags source_flags)
{
  if (source_flags & NO_FLAGS_SET) return;

  const size_t source_event_tag = lbs_solver_.GetSourceEventTag();
  Chi::log.LogEvent(source_event_tag, chi::ChiLog::EventType::EVENT_BEGIN);

  apply_fixed_src_       = (source_flags & APPLY_FIXED_SOURCES);
  apply_wgs_scatter_src_ = (source_flags & APPLY_WGS_SCATTER_SOURCES);
  apply_ags_scatter_src_ = (source_flags & APPLY_AGS_SCATTER_SOURCES);
  apply_wgs_fission_src_ = (source_flags & APPLY_WGS_FISSION_SOURCES);
  apply_ags_fission_src_ = (source_flags & APPLY_AGS_FISSION_SOURCES);
  suppress_wg_scatter_src_ = (source_flags & SUPPRESS_WG_SCATTER);

  //================================================== Get group setup
  gs_i_ = static_cast<size_t>(groupset.groups_.front().id_);
  gs_f_ = static_cast<size_t>(groupset.groups_.back().id_);

  first_grp_ = static_cast<size_t>(lbs_solver_.Groups().front().id_);
  last_grp_ = static_cast<size_t>(lbs_solver_.Groups().back().id_);

  default_zero_src_.assign(lbs_solver_.Groups().size(), 0.0);

  const auto& cell_transport_views = lbs_solver_.GetCellTransportViews();
  const auto& matid_to_src_map = lbs_solver_.GetMatID2IsoSrcMap();

  const size_t num_moments = lbs_solver_.NumMoments();
  const auto& ext_src_moments_local = lbs_solver_.ExtSrcMomentsLocal();

  const auto& m_to_ell_em_map =
    groupset.quadrature_->GetMomentToHarmonicsIndexMap();

  //================================================== Loop over local cells
  const auto& grid = lbs_solver_.Grid();
  // Apply all nodal sources
  for (const auto& cell : grid.local_cells)
  {
    auto& transport_view = cell_transport_views[cell.local_id_];
    cell_volume_ = transport_view.Volume();

    //==================== Obtain xs
    const auto& xs = transport_view.XS();

      std::shared_ptr<chi_physics::IsotropicMultiGrpSource> P0_src = nullptr;
    if (matid_to_src_map.count(cell.material_id_) > 0)
      P0_src = matid_to_src_map.at(cell.material_id_);

    const auto& S = xs.TransferMatrices();
    const auto& F = xs.ProductionMatrix();
    const auto& precursors = xs.Precursors();
    const auto& nu_delayed_sigma_f = xs.NuDelayedSigmaF();

    //======================================== Loop over nodes
    const int num_nodes = transport_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      //=================================== Loop over moments
      for (int m = 0; m < static_cast<int>(num_moments); ++m)
      {
        unsigned int ell = m_to_ell_em_map[m].ell;

        size_t uk_map = transport_view.MapDOF(i, m, 0); //unknown map

        const double* phi = &phi_local[uk_map];

        //==================== Declare moment src
        if (P0_src and ell == 0)
          fixed_src_moments_ = P0_src->source_value_g_.data();
        else
          fixed_src_moments_ = default_zero_src_.data();

        if (lbs_solver_.Options().use_src_moments)
          fixed_src_moments_ = &ext_src_moments_local[uk_map];

        //============================= Loop over groupset groups
        for (size_t g = gs_i_; g <= gs_f_; ++g)
        {
          g_ = g;

          double rhs = 0.0;

          //============================== Apply fixed sources
          if (apply_fixed_src_) rhs += this->AddSourceMoments();

          //============================== Apply scattering sources
          if (ell < S.size())
          {
            const auto& S_ell = S[ell];
            //==================== Add Across GroupSet Scattering (AGS)
            if (apply_ags_scatter_src_)
              for (const auto& [_, gp, sigma_sm] : S_ell.Row(g))
                if (gp < gs_i_ or gp > gs_f_)
                  rhs += sigma_sm * phi[gp];

            //==================== Add Within GroupSet Scattering (WGS)
            if (apply_wgs_scatter_src_)
              for (const auto& [_, gp, sigma_sm] : S_ell.Row(g))
                if (gp >= gs_i_ and gp <= gs_f_)
                {
                  if (suppress_wg_scatter_src_ and g_ == gp) continue;
                  rhs += sigma_sm * phi[gp];
                }
          }

          //============================== Apply fission sources
          const bool fission_avail = ell == 0 and xs.IsFissionable();

          if (fission_avail)
          {
            const auto& F_g = F[g];
            if (apply_ags_fission_src_)
              for (size_t gp = first_grp_; gp <= last_grp_; ++gp)
                if (gp < gs_i_ or gp > gs_f_)
                  rhs += F_g[gp] * phi[gp];

            if (apply_wgs_fission_src_)
              for (size_t gp = gs_i_; gp <= gs_f_; ++gp)
                rhs += F_g[gp] * phi[gp];

            if (lbs_solver_.Options().use_precursors)
              rhs += this->AddDelayedFission(
                precursors, nu_delayed_sigma_f, &phi_local[uk_map]);
          }

          //============================== Add to destination vector
          destination_q[uk_map + g] += rhs;

        }//for g
      }//for m
    }//for dof i
  }//for cell

  AddAdditionalSources(groupset, destination_q, phi_local, source_flags);

  Chi::log.LogEvent(source_event_tag, chi::ChiLog::EventType::EVENT_END);
}

double SourceFunction::AddSourceMoments() const
{
  return fixed_src_moments_[g_];
}


//###################################################################
/**Adds delayed particle precursor sources.*/
double SourceFunction::
  AddDelayedFission(const PrecursorList &precursors,
                    const std::vector<double> &nu_delayed_sigma_f,
                    const double *phi) const
{
  double value = 0.0;
  if (apply_ags_fission_src_)
    for (size_t gp = first_grp_; gp <= last_grp_; ++gp)
      if (gp < gs_i_ or gp > gs_f_)
        for (const auto& precursor : precursors)
          value += precursor.emission_spectrum[g_] *
                   precursor.fractional_yield *
                   nu_delayed_sigma_f[gp] * phi[gp];

  if (apply_wgs_fission_src_)
    for (size_t gp = gs_i_; gp <= gs_f_; ++gp)
      for (const auto& precursor : precursors)
        value += precursor.emission_spectrum[g_] *
                 precursor.fractional_yield *
                 nu_delayed_sigma_f[gp] * phi[gp];

  return value;
}


//###################################################################
/**Adds point sources to the source moments.*/
void SourceFunction::
  AddPointSources(LBSGroupset &groupset,
                  std::vector<double> &destination_q,
                  const std::vector<double>&,
                  SourceFlags source_flags)
{
  const bool apply_fixed_src       = (source_flags & APPLY_FIXED_SOURCES);

  const auto& cell_transport_views = lbs_solver_.GetCellTransportViews();

  const auto gs_i = static_cast<size_t>(groupset.groups_.front().id_);
  const auto gs_f = static_cast<size_t>(groupset.groups_.back().id_);

  //================================================== Apply point sources
  if (not lbs_solver_.Options().use_src_moments and apply_fixed_src)
    for (const auto& point_source : lbs_solver_.PointSources())
    {
      const auto& info_list = point_source.ContainingCellsInfo();
      for (const auto& info : info_list)
      {
        auto& transport_view = cell_transport_views[info.cell_local_id];

        const auto& strength = point_source.Strength();
        const auto& node_weights = info.node_weights;
        const double vol_w = info.volume_weight;

        const int num_nodes = transport_view.NumNodes();
        for (int i = 0; i < num_nodes; ++i)
        {
          const size_t uk_map = transport_view.MapDOF(i, /*moment=*/0, /*grp=*/0);
          for (size_t g = gs_i; g <= gs_f; ++g)
            destination_q[uk_map + g] += strength[g] * node_weights[i] * vol_w;
        }//for node i
      }//for cell
    }//for point source
}

}//namespace lbs