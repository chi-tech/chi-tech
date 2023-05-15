#include "lbs_sweepchunk.h"

namespace lbs
{
// ##################################################################
/**Assembles the volumetric gradient term.*/
void LBSSweepChunk::VolumetricGradientTerm()
{
  const auto& G = *G_;

  for (int i = 0; i < cell_num_nodes_; ++i)
    for (int j = 0; j < cell_num_nodes_; ++j)
      Amat_[i][j] = omega_.Dot(G[i][j]);
}

// ##################################################################
/**Performs the integral over the surface of a face.*/
void LBSSweepChunk::UpwindSurfaceIntegrals(size_t f)
{
  const auto& M_surf_f = (*M_surf_)[f];
  const double mu = face_mu_values_[f];
  const size_t num_face_nodes = cell_mapping_->NumFaceNodes(f);
  for (int fi = 0; fi < num_face_nodes; ++fi)
  {
    const int i = cell_mapping_->MapFaceNode(f, fi);
    for (int fj = 0; fj < num_face_nodes; ++fj)
    {
      const int j = cell_mapping_->MapFaceNode(f, fj);

      const double* psi = sweep_surface_status_info_.GetUpwindPsi(fj);

      const double mu_Nij = -mu * M_surf_f[i][j];
      Amat_[i][j] += mu_Nij;
      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
        b_[gsg][i] += psi[gsg] * mu_Nij;
    } // for face node j
  }   // for face node i
}

// ##################################################################
/**Assembles angular sources and applies the mass matrix terms.*/
void LBSSweepChunk::AssembleMassTerms(const int g,
                                      const int gsg,
                                      const double sigma_tg)
{
  const auto& M = *M_;
  const auto& m2d_op = groupset_.quadrature_->GetMomentToDiscreteOperator();

  // ============================= Contribute source moments
  // q = M_n^T * q_moms
  for (int i = 0; i < cell_num_nodes_; ++i)
  {
    double temp_src = 0.0;
    for (int m = 0; m < num_moments_; ++m)
    {
      const size_t ir = cell_transport_view_->MapDOF(i, m, g);
      temp_src += m2d_op[m][direction_num_] * q_moments_[ir];
    } // for m
    source_[i] = temp_src;
  } // for i

  // ============================= Mass Matrix and Source
  // Atemp  = Amat + sigma_tgr * M
  // b     += M * q
  for (int i = 0; i < cell_num_nodes_; ++i)
  {
    double temp = 0.0;
    for (int j = 0; j < cell_num_nodes_; ++j)
    {
      const double Mij = M[i][j];
      Atemp_[i][j] = Amat_[i][j] + Mij * sigma_tg;
      temp += Mij * source_[j];
    } // for j
    b_[gsg][i] += temp;
  } // for i
}

// ##################################################################
/**Adds a single direction's contribution to the moment integrals.*/
void LBSSweepChunk::ContributeDirectionToFluxMomentIntgls()
{
  const auto& d2m_op = groupset_.quadrature_->GetDiscreteToMomentOperator();

  auto& output_phi = GetDestinationPhi();

  for (int m = 0; m < num_moments_; ++m)
  {
    const double wn_d2m = d2m_op[m][direction_num_];
    for (int i = 0; i < cell_num_nodes_; ++i)
    {
      const size_t ir = cell_transport_view_->MapDOF(i, m, gs_gi_);
      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
        output_phi[ir + gsg] += wn_d2m * b_[gsg][i];
    }
  }
}

// ##################################################################
/**Updates angular fluxes.*/
void LBSSweepChunk::UpdateAngularFluxes()
{
  typedef const int64_t cint64_t;

  auto& output_psi = GetDestinationPsi();

  for (int i = 0; i < cell_num_nodes_; ++i)
  {
    cint64_t imap = grid_fe_view_.MapDOFLocal(
      *cell_, i, groupset_.psi_uk_man_, direction_num_, gs_ss_begin_);
    for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
      output_psi[imap + gsg] = b_[gsg][i];
  } // for i
}

// ##################################################################
/**Operations when outgoing fluxes are handled including passing
 * face angular fluxes downstream and computing
 * balance parameters (i.e. outflow)
 * */
void LBSSweepChunk::OutgoingSurfaceOperations(const size_t f)
{
  const auto& IntF_shapeI = (*IntS_shapeI_)[f];
  const double mu = face_mu_values_[f];
  const double wt = direction_qweight_;

  const auto& sss_info = sweep_surface_status_info_;

  const size_t num_face_nodes = cell_mapping_->NumFaceNodes(f);
  for (int fi = 0; fi < num_face_nodes; ++fi)
  {
    const int i = cell_mapping_->MapFaceNode(f, fi);

    double* psi = sweep_surface_status_info_.GetDownwindPsi(fi);

    if (not sss_info.on_boundary or sss_info.is_reflecting_bndry_)
      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
        psi[gsg] = b_[gsg][i];
    if (sss_info.on_boundary and not sss_info.is_reflecting_bndry_)
      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
        cell_transport_view_->AddOutflow(gs_gi_ + gsg,
                                         wt * mu * b_[gsg][i] * IntF_shapeI[i]);
  } // for fi
}

} // namespace lbs