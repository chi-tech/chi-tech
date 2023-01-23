#ifndef CHITECH_CHI_FFINTER_POINT_H
#define CHITECH_CHI_FFINTER_POINT_H

#include "../chi_ffinterpolation.h"
#include "ChiMesh/chi_mesh.h"

namespace chi_mesh
{
//###################################################################
/** A line based interpolation function.*/
  class FieldFunctionInterpolationPoint :
    public FieldFunctionInterpolation
  {
  public:
    chi_mesh::Vector3 m_point_of_interest;

  protected:
    bool m_locally_owned = false;
    uint64_t m_owning_cell_gid = 0;
    double m_point_value = 0.0;

  public:
    FieldFunctionInterpolationPoint() :
      FieldFunctionInterpolation(ff_interpolation::Type::POINT)
    {  }

    void Initialize() override;
    void Execute() override;
    double GetPointValue() const;

  public:
    std::string GetDefaultFileBaseName() const override
    {return "";}
    void ExportPython(std::string base_name) override {};
  };
}//namespace chi_mesh

#endif //CHITECH_CHI_FFINTER_POINT_H
