#ifndef CHI_FLUDS_H
#define CHI_FLUDS_H

#include <vector>
#include <set>
#include <cstddef>
#include <cstdint>

#include "FLUDSCommonData.h"

namespace chi_mesh
{
class GridFaceHistogram;

} // namespace chi_mesh

namespace chi_mesh::sweep_management
{

class SPDS;

class FLUDS
{
public:
  explicit FLUDS(size_t num_groups, size_t num_angles, const SPDS& spds)
    : num_groups_(num_groups), num_angles_(num_angles), spds_(spds){};

  const SPDS& GetSPDS() const { return spds_; }

  virtual double* OutgoingPsi(int cell_so_index,
                              int outb_face_counter,
                              int face_dof,
                              int n) = 0;
  virtual double* UpwindPsi(
    int cell_so_index, int inc_face_counter, int face_dof, int g, int n) = 0;

  virtual double*
  NLOutgoingPsi(int nonl_outb_face_counter, int face_dof, int n) = 0;

  virtual double*
  NLUpwindPsi(int nonl_inc_face_counter, int face_dof, int g, int n) = 0;

  virtual void ClearLocalAndReceivePsi() {}
  virtual void ClearSendPsi() {}
  virtual void AllocateInternalLocalPsi(size_t num_grps, size_t num_angles) {}
  virtual void
  AllocateOutgoingPsi(size_t num_grps, size_t num_angles, size_t num_loc_sucs)
  {
  }
  virtual void
  AllocateIncomingPsi(size_t num_grps, size_t num_angles, size_t num_loc_deps)
  {
  }
  virtual void AllocateDelayedLocalPsi(size_t num_grps, size_t num_angles) {}
  virtual void AllocatePrelocIOutgoingPsi(size_t num_grps,
                                          size_t num_angles,
                                          size_t num_loc_deps)
  {
  }
  virtual void AllocateDelayedPrelocIOutgoingPsi(size_t num_grps,
                                                 size_t num_angles,
                                                 size_t num_loc_deps)
  {
  }

  virtual std::vector<double>& DelayedLocalPsi() = 0;
  virtual const std::vector<double>& DelayedLocalPsi() const = 0;
  virtual std::vector<double>& DelayedLocalPsiOld() = 0;
  virtual const std::vector<double>& DelayedLocalPsiOld() const = 0;

  virtual std::vector<std::vector<double>>& DeplocIOutgoingPsi() = 0;
  virtual const std::vector<std::vector<double>>&
  DeplocIOutgoingPsi() const = 0;

  virtual std::vector<std::vector<double>>& PrelocIOutgoingPsi() = 0;
  virtual const std::vector<std::vector<double>>&
  PrelocIOutgoingPsi() const = 0;

  virtual std::vector<std::vector<double>>& DelayedPrelocIOutgoingPsi() = 0;
  virtual const std::vector<std::vector<double>>&
  DelayedPrelocIOutgoingPsi() const = 0;
  virtual std::vector<std::vector<double>>& DelayedPrelocIOutgoingPsiOld() = 0;
  virtual const std::vector<std::vector<double>>&
  DelayedPrelocIOutgoingPsiOld() const = 0;

  virtual ~FLUDS() = default;

protected:
  const size_t num_groups_;
  const size_t num_angles_;
  const SPDS& spds_;
};

} // namespace chi_mesh::sweep_management

#endif // CHI_FLUDS_H