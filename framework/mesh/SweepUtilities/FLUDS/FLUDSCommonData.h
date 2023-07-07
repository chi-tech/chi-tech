#ifndef CHITECH_FLUDSCOMMONDATA_H
#define CHITECH_FLUDSCOMMONDATA_H

#include <vector>

namespace chi_mesh::sweep_management
{
class SPDS;
struct FaceNodalMapping
{
  const int associated_face;
  const std::vector<short> node_mapping;

  FaceNodalMapping(int in_ass_face, std::vector<short>& in_node_mapping)
    : associated_face(in_ass_face), node_mapping(in_node_mapping)
  {
  }
};
typedef std::vector<FaceNodalMapping> CellFaceNodalMapping;

class FLUDSCommonData
{
public:
  explicit FLUDSCommonData(
    const SPDS& spds,
    const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

  virtual ~FLUDSCommonData() = default;

  const SPDS& GetSPDS() const;

protected:
  const SPDS& spds_;
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings_;
};

} // namespace chi_mesh::sweep_management

#endif // CHITECH_FLUDSCOMMONDATA_H
