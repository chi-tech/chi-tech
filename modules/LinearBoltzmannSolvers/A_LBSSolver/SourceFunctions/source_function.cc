#include "source_function.h"

#include "A_LBSSolver/lbs_solver.h"

#include "source_context.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace lbs
{

//###################################################################
/**Constructor.*/
SourceFunction::SourceFunction(const LBSSolver &lbs_solver,
                               std::shared_ptr<SourceContext>& context) :
  lbs_solver_(lbs_solver),
  grid_(lbs_solver_.Grid()),
  context_(context)
{}


//###################################################################
/**Sets the source moments for the groups in the current group set.
 *
 * \param groupset The groupset the under consideration.
 * \param destination_q A vector to contribute the source to.
 * \param phi The primary STL vector to operate off.
 * \param source_flags Flags for adding specific terms into the
 *        destination vector. Available flags are for applying
 *        the material source, across/within-group scattering,
 *        and across/within-groups fission.
 *
 * */
void SourceFunction::operator()(LBSGroupset &groupset,
                                std::vector<double> &destination_q,
                                const std::vector<double> &phi,
                                SourceFlags source_flags)
{
  const size_t source_event_tag = lbs_solver_.GetSourceEventTag();
  chi::log.LogEvent(source_event_tag, chi_objects::ChiLog::EventType::EVENT_BEGIN);

  const bool apply_fixed_src       = (source_flags & APPLY_FIXED_SOURCES);
  const bool apply_wgs_scatter_src = (source_flags & APPLY_WGS_SCATTER_SOURCES);
  const bool apply_ags_scatter_src = (source_flags & APPLY_AGS_SCATTER_SOURCES);
  const bool apply_wgs_fission_src = (source_flags & APPLY_WGS_FISSION_SOURCES);
  const bool apply_ags_fission_src = (source_flags & APPLY_AGS_FISSION_SOURCES);
  const bool suppress_wg_scatter_src = (source_flags & SUPPRESS_WG_SCATTER);

  context_->SetFlags(apply_wgs_scatter_src,
                     apply_ags_scatter_src,
                     apply_wgs_fission_src,
                     apply_ags_fission_src,
                     suppress_wg_scatter_src);

  //================================================== Get group setup
  const auto gs_i = static_cast<size_t>(groupset.groups_.front().id_);
  const auto gs_f = static_cast<size_t>(groupset.groups_.back().id_);

  context_->SetGroupsetBounds(gs_i, gs_f);

  const auto first_grp = static_cast<size_t>(lbs_solver_.Groups().front().id_);
  const auto last_grp = static_cast<size_t>(lbs_solver_.Groups().back().id_);

  context_->SetGroupBounds(first_grp, last_grp);

  std::vector<double> default_zero_src(lbs_solver_.Groups().size(), 0.0);

  const auto& cell_transport_views = lbs_solver_.GetCellTransportViews();
  const auto& matid_to_src_map = lbs_solver_.GetMatID2IsoSrcMap();

  const size_t num_moments = lbs_solver_.NumMoments();
  const auto& ext_src_moments_local = lbs_solver_.ExtSrcMomentsLocal();

  const auto& m_to_ell_em_map =
    groupset.quadrature_->GetMomentToHarmonicsIndexMap();

  //================================================== Loop over local cells
  // Apply all nodal sources
  for (const auto& cell : grid_.local_cells)
  {
    auto& transport_view = cell_transport_views[cell.local_id_];
    context_->SetCellVolume(transport_view.Volume());

    //==================== Obtain xs
    const auto& xs = transport_view.XS();
    auto P0_src = matid_to_src_map.at(cell.material_id_);

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

        //==================== Declare moment src
        const double* src;
        if (P0_src and ell == 0)
          src = P0_src->source_value_g_.data();
        else
          src = default_zero_src.data();

        if (lbs_solver_.Options().use_src_moments)
          src = &ext_src_moments_local[uk_map];

        context_->SetFixedSrcMomentsData(src);

        //============================= Loop over groupset groups
        for (size_t g = gs_i; g <= gs_f; ++g)
        {
          context_->SetGroupIndex(g);

          double rhs = 0.0;

          //============================== Apply fixed sources
          if (apply_fixed_src) rhs += context_->AddSourceMoments();

          //============================== Apply scattering sources
          if (ell < S.size()) rhs += context_->AddScattering(S[ell], &phi[uk_map]);

          //============================== Apply fission sources
          const bool fission_avail = ell == 0 and xs.IsFissionable();

          if (fission_avail)
          {
            rhs += context_->AddPromptFission(F[g], &phi[uk_map]);
            if (lbs_solver_.Options().use_precursors)
              rhs += context_->AddDelayedFission(
                precursors, nu_delayed_sigma_f, &phi[uk_map]);
          }

          //============================== Add to destination vector
          destination_q[uk_map + g] += rhs;

        }//for g
      }//for m
    }//for dof i
  }//for cell

  AddAdditionalSources(groupset, destination_q, phi, source_flags);

  chi::log.LogEvent(source_event_tag, chi_objects::ChiLog::EventType::EVENT_END);
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