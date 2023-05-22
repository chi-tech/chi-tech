#ifndef CELL_MAPPING_FE_PWL_BASE_H
#define CELL_MAPPING_FE_PWL_BASE_H

#include <utility>

#include "ChiMesh/chi_mesh.h"

#include "ChiMath/SpatialDiscretization/CellMappings/cell_mapping_base.h"

//###################################################################
namespace chi_math
{
  /** Base class for all cell FE views.*/
  class CellMappingFE_PWL : public CellMapping
  {
  protected:
    const std::vector<chi_mesh::Vector3> node_locations_;

  public:
    typedef std::vector<double> VecDbl;
    typedef std::vector<chi_mesh::Vector3> VecVec3;

  public:
    /** Constructor. */
    explicit CellMappingFE_PWL(const chi_mesh::MeshContinuum& ref_grid,
                               const chi_mesh::Cell& cell,
                               size_t num_nodes,
                               std::vector<chi_mesh::Vector3> vertices,
                               std::vector<std::vector<int>> face_node_mappings) :
        CellMapping(ref_grid, cell, num_nodes,
                    std::move(face_node_mappings),
                    &CellMapping::ComputeCellVolumeAndAreas),
        node_locations_(std::move(vertices))
    {}

    //02 ShapeFuncs
    std::vector<chi_mesh::Vector3> GetNodeLocations() const override;

  protected:
    /** Spatial weight function. See also ComputeWeightedUnitIntegrals. */
    virtual double SpatialWeightFunction(const chi_mesh::Vector3& pt) const
    { return 0.0; }

    /** Compute spatially-weighted unit integrals from quadrature point data.
     *  Spatial weighting is given by SpatialWeightFunction evaluated at
     *  the real quadrature points. */
    void ComputeWeightedUnitIntegrals(
        finite_element::UnitIntegralData& ui_data) const;

  protected:
    static std::vector<chi_mesh::Vector3>
    GetVertexLocations(const chi_mesh::MeshContinuum& grid,
                       const chi_mesh::Cell& cell);

    static std::vector<std::vector<int>>
    MakeFaceNodeMapping(const chi_mesh::Cell& cell);
  };
}


#endif